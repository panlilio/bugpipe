function segmentAndTrackMP_spmd(D,doSegment,doTrack)
%SEGMENTANDTRACKMP_SPMD. 
%2016-10-29. Segments and tracks cells for the dataset with directory
%structure D (obtained from either getDirsMP or getDirsMP_SSD), 
%using PARFOR and SPMD from the Parallel Computing Toolbox.

totTic = tic;

%Get the orientation of the device, store in tracking parameters
dataInfo = lookupDataMP(D.baseDir);
P = segmentationParamsMP(dataInfo.ORIENTATION);

%Retrieve png images (including any blanks, if they exist)
L = imListMP(D.dataDir);

%Create a diary to record processing times
diaryMP = fullfile(D.lineDir,'diary.txt');
if exist(diaryMP,'file')
    fID = fopen(diaryMP,'a');
else
    fID = fopen(diaryMP,'w');
end
fprintf(fID,'%s \t started segmentation. \r',datestr(clock));
fclose(fID);

%Variables to pass into parfor and spmd
dataDir = D.dataDir;
matDir = D.matDir;
mat2Dir = D.mat2Dir;
lineDir = D.lineDir;

IM = L.imList;
BL = L.blankList;
FT = L.imListFOVT;
TA = L.timeCreated;

%% SEGMENTATION (FOV AND TIME-INDEPENDENT)
segTic = tic;

if doSegment
    spmd(10)
        imList = codistributed(IM,codistributor1d(1));
        blankList = codistributed(BL,codistributor1d(1));
        fovT = codistributed(FT,codistributor1d(1));
        tAcq = codistributed(TA,codistributor1d(1));

        im = getLocalPart(imList);
        bl = getLocalPart(blankList);
        fT = getLocalPart(fovT);
        tA = getLocalPart(tAcq);

        for k = 1:numel(im)
            %Read in image and blank
            I = imread(fullfile(dataDir,im{k}));
            if ~isempty(bl{k})
                Iblank = imread(fullfile(dataDir,bl{k}));
            else
                Iblank = zeros(512);
            end

            I = double(I) - double(Iblank);
            I(I<0) = 0;

            %Segment
            [CELLS,CHANNELS] = segmentationMP(I,P);

            %Record actual time observed
            CELLS.tAcq = tA(k);

            %Save variables
            parsave1(fullfile(matDir,sprintf('fov%.2d_t%.4d',fT(k,1),fT(k,2))),CELLS,CHANNELS)
        end
    end
end

segToc = toc(segTic);

fID = fopen(diaryMP,'a');
fprintf(fID,'%s \t finished segmentation. \r \t Segmentation time: %f sec \r',datestr(clock),segToc);
fprintf(fID,'%s \t started tracking. \r',datestr(clock));
fclose(fID);

%% TRACK (FOV-INDEPENDENT)
fovRange = L.fovRange;
tRange = L.tRange;
deadEnd = P.deadEnd;

trackTic = tic;

if doTrack
    parfor f = 1:numel(L.fovRange)
        fov = fovRange(f);
        C0 = {};
        I0 = zeros(512);
        for t = tRange 
            %Load the current segmented data
            mat1 = fullfile(matDir,sprintf('fov%.2d_t%.4d.mat',fov,t)); 
            if ~exist(mat1,'file'), 
                fprintf('\tCONTINUE: No MAT found. fov = %d, t = %d \r',fov,t) 
                continue
            end
            C1 = load(mat1);

            if isempty(C1.CELLS.boundary) || isempty(C1.CHANNELS)
                fprintf('\tCONTINUE: No cells segmented. fov = %d, t = %d \r',fov,t) 
                continue
            end

            %Read in the previous and current image
            [imname,blankname] = L.getfilenamesFOVT(fov,t,1); %#ok<PFBNS>
            I1 = imread(imname);
            if isempty(blankname)
                Iblank = zeros(size(I1));
            else
                Iblank = imread(blankname);
            end

            I1 = double(I1) - double(Iblank);
            I1(I1<0) = 0;

            %Track cells
            [C2,~] = cellTrackerMP(I1,I0,C1,C0,t,fov,lineDir,deadEnd);

            %Save variables
            CELLS = C2.CELLS;
            CHANNELS = C2.CHANNELS;
            parsave1(fullfile(mat2Dir,sprintf('fov%.2d_t%.4d',fov,t)),CELLS,CHANNELS)

            %Update for the next iteration but skip empty frames from tracking
            %consideration.
            if ~isempty(CELLS.boundary)
                C0 = C2;
                I0 = I1;
            end
        end
    end
end

trackToc = toc(trackTic);
totToc = toc(totTic);

fID = fopen(diaryMP,'a');
fprintf(fID,'%s \t finished tracking. \r \t Tracking time: %f sec \r',datestr(clock),trackToc);
fprintf(fID,'nFOV \t %d \r',numel(L.fovRange));
fprintf(fID,'nTimepoints \t %d \r',numel(L.tRange));
fprintf(fID,'nImages \t %d \r',numel(L.imList));
fprintf(fID,'Segmentation speed \t %f fps \r',numel(L.imList)/segToc);
fprintf(fID,'Tracking speed \t %f fps \r',numel(L.imList)/trackToc);
fprintf(fID,'Total speed \t %f fps \r',numel(L.imList)/totToc);
fclose(fID);

end