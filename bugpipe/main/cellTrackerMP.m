function [C2,STATS] = cellTrackerMP(I1,I0,C1,C0,t,fov,dirname,deadEnd)
%[C2,STATS] = CELLTRACKERMP(I1,C1,C0,t,fov,dirname,deadEnd). 
%Tracks cells and cell lines from the mother machine. Obtains morphology and 
%fluorescence properties of each cell. Saves and reads in cell lines as MATs. 
%Input I1 and C1 pertain to the frame of interest, C0 refers to the previous 
%frame. Matrices C0 and C1 are the loaded structured arrays with fields CELLS 
%and CHANNELS, the output from segmentationMP. Cell line numbering convention 
%starts at 000001 for each field of view, i.e. need field of view number to 
%properly save and load data. Output C2 is the updated structure containing 
%CELLS and CHANNELS from C1; STATS keeps track of the cell tracking errors.

%Tracking and labelling parameters
strFormatLabel = '%.2d_%.4d_%.9d_0'; %fov,channel,cellno. Base cell names. Division markers are appended.
strFormatLine = '%.2d_%.4d_INFO'; %fov,channel. Base name for entire cell line.
strFormatFOVDir = 'fov%.2d'; %format for naming the field of view directory.
strFormatLineDir = 'line%.4d'; %format for naming the line (channel) directory.

%Remove cells touching the border
C1 = rmBorderCellsMP(C1);

%Skip tracking if all objects were touching the border
if isempty(C1.CELLS.boundary) 
    C2 = C1;
    STATS = {};
    return
end

%Initialize variables
STATS = {};

CELLS1 = C1.CELLS;
CHANNELS1 = C1.CHANNELS;
nCells1 = numel(CELLS1.area);
nChannels1 = numel(CHANNELS1.threshold);
CELLS1.rank = zeros(1,nCells1);
CELLS1.label = cell(nCells1,1);
CELLS1.lineID = zeros(1,nCells1);

if isempty(C0) 
    %No previous image to work off of--first image in the time series. Label all cells naively. 
    CHANNELS1.label = sort(unique(CELLS1.channel));
    channelLabel = CHANNELS1.label;
    
    for k = CHANNELS1.label 
        %Directory to hold MAT cell data    
        linedirname = getLineDirnameMP(k);
        
        %Label cells naively
        [cellLine,cellLineFrames,cellLineFamilyMap] = naiveLabelMP(channelLabel(k),linedirname);  %#ok<ASGLU>
        
        %Save the cell line variables to access later
        matname = sprintf(strFormatLine,fov,k);
        save(fullfile(linedirname,matname),'cellLine','cellLineFrames','cellLineFamilyMap','-v6')
    end
    
    C2 = {};
    C2.CHANNELS = CHANNELS1;
    C2.CELLS = CELLS1;
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%IF THERE IS A PREVIOUS FRAME TO USE IN TRACKING

%Variables for the previous frame
CELLS0 = C0.CELLS;
CHANNELS0 = C0.CHANNELS;
nChannels0 = numel(CHANNELS0.threshold);

%Find nearest channel in the previous frame, will also update CHANNELS1.label and
%CELLS1.channel
[channelNN,channelLabel] = nearestChannelMP;

%Perform rank comparison of cells
for k = 1:nChannels1
    %Cell line MAT variable
    matname = sprintf(strFormatLine,fov,channelLabel(k));
    matname = strcat(matname,'.mat');
    
    linedirname = getLineDirnameMP(channelLabel(k));
        
    if isfinite(channelNN(k)) && exist(fullfile(linedirname,matname),'file')
        %Label using the nearest channel in the previous image
        [cID,LID,parentID] = rankLabelMP(channelLabel(k));
        
        %Store cell properties in relevant indices of cellLine
        M = load(fullfile(linedirname,matname)); 
        cellLine = M.cellLine; 
        cellLineFrames = M.cellLineFrames; 
        cellLineFamilyMap = M.cellLineFamilyMap; 
        
        lc = numel(M.cellLine)+1;
        for c = 1:numel(cID)
            cellID = cID(c);
            P = getCellPropsMP(CELLS1.boundary{cellID},I1,CELLS1.tAcq);           
            if LID(c)==0
                %New entry in cellLine required
                cellData = cellDataMP;
                cellData.frames = [t t];
                cellData.label = CELLS1.label{cellID};
                cellData.cellID = [cellData.cellID cellID];
                cellData.lineID = lc;
                cellData.parent = parentID(c);
                
                cellData = cellData.saveProps(P);
                
                %Save cell data
                save(fullfile(linedirname,cellData.label),'cellData','-v6')
                
                %Update                
                cellLine{lc} = cellData.label;
                cellLineFrames = [cellLineFrames; [t t]]; %#ok<AGROW>
                cellLineFamilyMap = [cellLineFamilyMap parentID(c)]; %#ok<AGROW>
                
                CELLS1.lineID(cellID) = lc;
                lc = lc + 1;
                
                if parentID(c)~=0
                    %Division occurred
                    %Update parent
                    parentLabel = cellLine{parentID(c)}; 
                    cellData = load(fullfile(linedirname,strcat(parentLabel,'.mat')));
                    cellData = cellData.cellData;
                    cellData.divisionTime = diff(cellData.frames)+1;
                    
                    %Save updated parent
                    save(fullfile(linedirname,cellData.label),'cellData','-v6')
                end
            else
                %Update existing element in cellLine
                cellLabel = CELLS1.label{cellID};
                cellData = load(fullfile(linedirname,strcat(cellLabel,'.mat')));
                cellData = cellData.cellData;
                cellData.frames(2) = t;
                cellData.cellID = [cellData.cellID cellID];
                cellData = cellData.saveProps(P);
                
                %Save updated cell
                save(fullfile(linedirname,strcat(cellData.label,'.mat')),'cellData','-v6')
                                
                cellLine{LID(c)} = cellData.label; 
                cellLineFrames(LID(c),2) = t; 
                
                CELLS1.lineID(cellID) = LID(c);
            end                
        end
    else
        %Label cells naively 
        [cellLine,cellLineFrames,cellLineFamilyMap] = naiveLabelMP(channelLabel(k),linedirname); %#ok<ASGLU>
    end
    
    %Save the cell line variables to access later
    save(fullfile(linedirname,matname),'cellLine','cellLineFrames','cellLineFamilyMap','-v6')
end

%Output is updated structure containing CELLS1 and CHANNELS1
C2 = {};
C2.CELLS = CELLS1;
C2.CHANNELS = CHANNELS1;




%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ADDITIONAL FUNCTIONS CALLED BY MAIN BODY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function linedirname = getLineDirnameMP(chLabel)
%Construct the directory name according to channel and field of view
    linedirname = fullfile(dirname,sprintf(strFormatFOVDir,fov),sprintf(strFormatLineDir,chLabel));
    if ~exist(linedirname,'dir')
        mkdir(linedirname)
    end
end

%%
function R = rankCellsMP(xCOM)
%R = RANKCELLSMP(xCOM,deadend). Ranks cells according to x-centroid.

switch deadEnd
    case 'left'
        sortMode = 'ascend';
    case 'right'
        sortMode = 'descend';
end
[~,sid] = sort(xCOM,1,sortMode);
[~,ssid] = sort(sid);
R = ssid;

end %function rankCELLSMP


%%
function [cellLine,cellLineFrames,cellLineFamilyMap] = naiveLabelMP(chLabel,linedirname)
%Outputs cell line variables from nascent cell lines, nascent either
%because it is the first image of a stack or a channel was not tracked
%well. CELLS1 is updated with the field 

%Get index of cells of interest
cID1 = CELLS1.channel==chLabel;

%Rank cells
R1 = rankCellsMP(CELLS1.centroid(cID1,1));
CELLS1.rank(cID1) = R1;

%Initialize cell line array to save cell properties
cellLine = {};
cellLineFrames = []; %First and last frame for each cell in cellLine
cellLineFamilyMap = []; %Parent index of each cell in cellLine
jc = 1;
for j = find(cID1)
    %Label cells by field of view, channel, and rank
    CELLS1.label{j} = sprintf(strFormatLabel,fov,chLabel,R1(jc));
    CELLS1.lineID(j) = jc;
    
    %Get cell properties and store them
    P = getCellPropsMP(CELLS1.boundary{j},I1,CELLS1.tAcq);
    cellData = cellDataMP;
    cellData = cellData.saveProps(P);
    cellData.frames = [t t];
    cellData.label = CELLS1.label{j};
    cellData.lineID = jc;
    cellData.cellID = [cellData.cellID j];
       
    %Save MAT
    save(fullfile(linedirname,cellData.label),'cellData','-v6')
    
    %Update cellLine variables
    cellLine{jc} = cellData.label; %#ok<AGROW>
    cellLineFrames = [cellLineFrames; [t t]]; %#ok<AGROW> %First and last frame of each cell
    cellLineFamilyMap = [cellLineFamilyMap 0]; %#ok<AGROW>
    jc = jc + 1;
end

end %function naiveLabelMP


%%
function [channelNN,channelLabel] = nearestChannelMP
%Outputs the index of the nearest neighbour channel the previous frame 
%(index consistent with CHANNELS0.boundingBox order), and the channel 
%label that is consistent with between frames. Will update the labels in
%both CHANNELS1 and CELLS1.

%Remove drift by translating each frame so that the centre of mass of 
%the fluorescent image is at the centre
frameCOM0 = sum(I0,2);
frameCOM0 = frameCOM0./sum(frameCOM0);
frameCOM0 = sum(frameCOM0.*(1:size(I0,1))');

frameCOM1 = sum(I1,2);
frameCOM1 = frameCOM1./sum(frameCOM1);
frameCOM1 = sum(frameCOM1.*(1:size(I1,1))');

bb0 = mean(CHANNELS0.boundingBox(:,[1 2]),2);
bb1 = mean(CHANNELS1.boundingBox(:,[1 2]),2);
dChannel = (bb0-frameCOM0)*ones(1,nChannels1) - ones(nChannels0,1)*(bb1'-frameCOM1);
dChannel = dChannel.^2;

%Nearest channel in the previous frame to the current channels
mindChannel1 = min(dChannel,[],1);
nnChannelPrev_ = ones(nChannels0,1)*mindChannel1;
[nnChannelPrev1,channelID1] = find(dChannel==nnChannelPrev_);
channels1 = 1:nChannels1;

%Nearest channels in the current frame to the previous channels
mindChannel0 = min(dChannel,[],2);
nnChannelPrev_ = mindChannel0*ones(1,nChannels1);
[nnChannelPrev0,channelID0] = find(dChannel==nnChannelPrev_);

%Check for mutually nearest neighbours
NN = zeros(0,2);
for nc = 1:nChannels0
    c01 = channelID0(nnChannelPrev0==nc); %channel in current frame whose nearest neighbour is the ncth channel
    c10 = channelID1(nnChannelPrev1==nc); %nearest neighbour of the ncth channel in previous frame
    if any(ismember(c01,c10))
        f = ismember(c01,c10);
        f = find(f,1,'first');
        NN = [NN; nc c01(f)]; %#ok<AGROW>
    end
end

%In case some channels were omitted from the mutual nearest neighbour pairs
for nc = find(not(ismember(channels1,NN(:,2))))
    NN = [NN; NaN nc]; %#ok<AGROW>
end
[~,sid] = sort(NN(:,2));
NN = NN(sid,:); %Order pairs by index in frame1
channelNN = NN(:,1); %Channel index in frame0 of the nearest channel

%Get channel labels (can be different from index, because index indicates
%order of appearance along image y-axis, while label is maintained
%frame-to-frame, independent of e.g. drift in fields of view).
channelLabel = 0*channels1;
channelsOld = CELLS1.channel;
for nc = 1:nChannels1
    if isnan(NN(nc,1))
        channelLabel(nc) = randi([100 9999],1); %This will determine the MAT file label
    else
        channelLabel(nc) = CHANNELS0.label(channelNN(nc));
    end
    
    %Update label in CELLS1
    CELLS1.channel(channelsOld==nc) = channelLabel(nc);
end

%Update label in CHANNELS
CHANNELS1.label = channelLabel;

end %function nearestChannelMP



%%
function [cID1,LID,parentID] = rankLabelMP(chLabel)
%Tracks and labels cells in the logical index cID1 according to cells in
%the previous image, indexed by cID0. LID is the index in cellLine to
%store cell properties. LID is 0 if a new cell must be added to the line.
%parentID indicates the index in cellLine of the parent of a newly divided
%cell.

doLabel = 1;
resegmentationDone = 0;
divThresh = 0.7;

while doLabel 
    cID1 = find(CELLS1.channel==chLabel);    
    cID0 = find(CELLS0.channel==chLabel);
    
    %Get relevant properties
    R1 = rankCellsMP(CELLS1.centroid(cID1,1));
    CELLS1.rank(cID1) = R1;
    R0 = CELLS0.rank(cID0);
    
    [R1,sid1] = sort(R1);
    [R0,sid0] = sort(R0);
    
    A1 = CELLS1.area(cID1);
    A1 = A1(sid1);
    A0 = CELLS0.area(cID0);
    A0 = A0(sid0);
    
    cID1 = cID1(sid1);
    cID0 = cID0(sid0);
    
    %Initialize variables
    LID = zeros(1,numel(cID1));
    parentID = LID;
    r0 = 1;
    r1 = 1;
    
    %Step through ranks, checking area changes
    while r1 <= max(R1) 
        if any(R0==r0)
            %There is a cell in the previous frame to look at
            a1 = A1(r1);
            a0 = A0(r0);
            
            if a1 <= divThresh*a0
                %%%%%%%%%%%%%%%%%%%%%%DIVISION FLAG%%%%%%%%%%%%%%%%%%%%%%%%%%
                %Cell shrinks considerably.
                %Get area of next cell in line, if present
                if numel(A1)>=(r1+1)
                    a2 = A1(r1+1);
                else
                    a2 = [];
                end
                
                if isempty(a2) || (~isempty(a2) && a2<=divThresh*a0)
                    %At least one daughter cell is present
                    mID = cID0(r0); %mother index in CELLS0
                    mLabel = CELLS0.label{mID};
                    mLID = CELLS0.lineID(mID);
                    
                    CELLS1.label{cID1(r1)} = strcat(mLabel,'0');
                    parentID(r1) = mLID;
                    LID(r1) = 0;
                    
                    if ~isempty(a2)
                        %Second daughter cell is present as well
                        CELLS1.label{cID1(r1+1)} = strcat(mLabel,'1');
                        parentID(r1+1) = mLID;
                        LID(r1+1) = 0;
                        r1 = r1 + 1;
                    end
                elseif ~isempty(a2) && resegmentationDone<r1
                    %%%%%%%%%%%%%%%RESEGMENTATION FLAG 0%%%%%%%%%%%%%%%%%%%
                    %Next cell is too large to be a daughter cell: try
                    %segmenting it again.
                    tf = resegmentMP(cID1(r1+1));
                    resegmentationDone = r1;
                    if tf==1
                        %Resegmentation updated CELLS1 with new cells. Redo
                        %labelling of entire channel
                        break
                    else
                        %Labeling is unreliable. Need to initiate new cell
                        %labels for this cell and those remaining in the channel.
                        clck = clock; %Using the time to label ensures no repeated labels
                        labelN = [num2str(clck(4),'%.2d') num2str(clck(5),'%.2d') num2str(round(clck(6)*100),'%.5d')];
                        labelN = str2double(labelN) + 100;
                        for r = r1:max(R1)
                            CELLS1.label{cID1(r)} = sprintf(strFormatLabel,fov,chLabel,labelN);
                            parentID(r) = 0;
                            LID(r) = 0;
                            labelN = labelN + 1;
                        end
                        r1 = max(R1)+1;
                    end
                else
                    %CELL SHRINKS BUT SECOND DAUGHTER WAS SEGMENTED PROPERLY.
                    %Labeling is unreliable. Need to initiate new cell
                    %labels for this cell and those remaining in the channel.
                    clck = clock; %Using the time to label ensures no repeated labels
                    labelN = [num2str(clck(4),'%.2d') num2str(clck(5),'%.2d') num2str(round(clck(6)*100),'%.5d')];
                    labelN = str2double(labelN) + 100;
                    for r = r1:max(R1)
                        CELLS1.label{cID1(r)} = sprintf(strFormatLabel,fov,chLabel,labelN);
                        parentID(r) = 0;
                        LID(r) = 0;
                        labelN = labelN + 1;
                    end
                    r1 = max(R1)+1;
                end
                
                r0 = r0 + 1;
                r1 = r1 + 1;
                
            elseif a1 > 1.7*a0 && resegmentationDone<r1
                %%%%%%%%%%%%%%%%%%%%RESEGMENTATION FLAG%%%%%%%%%%%%%%%%%%%%%
                %Cell grows inordinately quickly
                tf = resegmentMP(cID1(r1)); 
                resegmentationDone = r1;
                if tf==1
                    %Resegmentation updated CELLS1 with new cells. Redo
                    %labelling of entire channel.
                    break
                else
                    %Labeling is unreliable. Need to initiate new cell
                    %labels for this cell and those remaining in the channel.
                    clck = clock; %Using the time to label ensures no repeated labels
                    labelN = [num2str(clck(4),'%.2d') num2str(clck(5),'%.2d') num2str(round(clck(6)*100),'%.5d')];
                    labelN = str2double(labelN) + 100;
                    for r = r1:max(R1)
                        CELLS1.label{cID1(r)} = sprintf(strFormatLabel,fov,chLabel,labelN);
                        parentID(r) = 0;
                        LID(r) = 0;
                        labelN = labelN + 1;
                    end
                    r1 = max(R1)+1;
                end
                
            else
                %%%%%%%%%%%%%%SIMPLE TRACKING: NO DIV OR RESEG%%%%%%%%%%%%%%
                mID = cID0(r0);
                mLabel = CELLS0.label{mID};
                mLID = CELLS0.lineID(mID);
                
                CELLS1.label{cID1(r1)} = mLabel;
                parentID(r1) = 0;
                LID(r1) = mLID;
                
                r0 = r0 + 1;
                r1 = r1 + 1;
            end
        else
            %There are too many cells in the current frame than can be
            %accounted for by the detected divisions. Might be a cell that
            %has arrived upstream in the device. Initiate a new label.
            %(This section identical to unreliable labelling from failed
            %resegmentation condition)
            clck = clock; %Using the time to label ensures no repeated labels
            labelN = [num2str(clck(4),'%.2d') num2str(clck(5),'%.2d') num2str(round(clck(6)*100),'%.5d')];
            labelN = str2double(labelN)+100; %Increase to avoid naive rank label repeats around midnight
            for r = r1:max(R1)
                CELLS1.label{cID1(r)} = sprintf(strFormatLabel,fov,chLabel,labelN);
                parentID(r) = 0;
                LID(r) = 0;
                labelN = labelN + 1;
            end
            r1 = max(R1)+1;
            
        end
    end %while r1<=max(R1)  
    
    if r1 > max(R1)
        doLabel = 0;
    end

end %while doLabel

end %function rankLabelMP


%%
function tf = resegmentMP(cellID)
%Resegments the cell in CELLS1 at index cellID, using erosions
%to isolate objects although the original boundary is retained.
%Structure CELLS is updated accordingly. Output tf is TRUE (1) 
%if resegmentation was successful

%Parameters
resizeScale = 3;
cellAreaMin = 30;

%TRY TO FIND ADDITIONAL CELLS BY BOUNDARY PROCESSING
CELL0 = getCellsMP(CELLS1.boundary{cellID}(:,[2 1]));
cellArea = zeros(numel(CELL0),1);
for nc = 1:numel(CELL0)
    A = getCellPropsMP(CELL0{nc});
    cellArea(nc) = A.area;
end

if all(cellArea>cellAreaMin)
    %UPDATE CELLS1
    %Clear all fields referring to the old boundary
    chLabel = CELLS1.channel(cellID);
    
    CELLS1.boundary(cellID) = [];
    CELLS1.area(cellID) = [];
    CELLS1.centroid(cellID,:) = [];
    CELLS1.channel(cellID) = [];
    CELLS1.rank(cellID) = []; %**nb. rank will be updated elsewhere (doLabel loop)
    
    if cellID<=numel(CELLS1.label)
        CELLS1.label(cellID) = [];
    end
    
    if cellID<=numel(CELLS1.lineID);
        CELLS1.lineID(cellID) = [];
    end
    
    nCellsLeft = numel(CELLS1.area);
    
    for nc = 1:numel(CELL0)
        A = getCellPropsMP(CELL0{nc});
        cellArea = A.area;
        if cellArea >= cellAreaMin
            CELLS1.boundary{nCellsLeft+1} = CELL0{nc};
            CELLS1.area(nCellsLeft+1) = cellArea;
            CELLS1.centroid(nCellsLeft+1,:) = mean(CELLS1.boundary{nCellsLeft+1});
            CELLS1.channel(nCellsLeft+1) = chLabel;
            nCellsLeft = nCellsLeft + 1;
        end
    end
    
    tf = 1;
    return
end

%IF BOUNDARY PROCESSING DOES NOT YIELD DISTINCT CELLS, TRY HARSHER
%THRESHOLD
%Get bounding box of image
xx = CELLS1.boundary{cellID}(:,1);
yy = CELLS1.boundary{cellID}(:,2);

xx = min(512,max(1,xx));
yy = min(512,max(1,yy));


bb = [floor(min(yy))-1 ceil(max(yy))+1 floor(min(xx))-1 ceil(max(xx))+1];
bb = min(512,max(1,bb));

Ibb = I1(bb(1):bb(2),bb(3):bb(4));
Ibb = medfilt2(Ibb);
Ireseg = imresize(Ibb,resizeScale);
Ireseg = mat2gray(Ireseg);

tOtsu = graythresh(Ireseg);
lowI = Ireseg(Ireseg(:)<tOtsu);
sig = sqrt(mean(lowI.^2)-mean(lowI).^2);
T = tOtsu + 2*sig;
BW = Ireseg>T;

%Reconstruct thresholded image using original boundary (on rescaled
%coordinates)
xx = (xx - bb(3))*resizeScale + 1; %Translate and scale to bounding box coords
yy = (yy - bb(1))*resizeScale + 1; %---------"----------

[xx,yy] = constructBoundaryMP(xx,yy);

cmyx = round(mean([yy xx],1));
J = false(size(BW));

JID = sub2ind(size(J),yy,xx);
J(JID) = true;
J = imfill(J,cmyx,4);
J = J~=0;

BW(not(J)) = 0;
L = bwlabel(BW);
nFound = max(L(:));

if nFound<=1
    %Resegmentation failed to produce separate cells.
    tf = 0;
    return
elseif nFound > 1
    %Using the original boundary, assign enclosed pixels to the nearest
    %object found during resegmentation.
        
    %Find the nearest newly found object to the points enclosed by the old
    %boundary
    [y,x] = find(J); %J must be on the same coordinate system ie bounding box as COM
    nPoints = numel(x);
    
    %Select random points to calculate distance from
    nRand = 64;
    dr = zeros(nPoints,nFound);
    for nf = 1:nFound
        [row,col] = find(L==nf);
        rID = randi(numel(row),nRand,1);
        xR = col(rID);
        yR = row(rID);
                
        dx = x*ones(1,nRand) - ones(nPoints,1)*(xR');
        dy = y*ones(1,nRand) - ones(nPoints,1)*(yR');
        dr(:,nf) = min(dx.^2+dy.^2,[],2) + 1e-10*(rand(size(dx,1),1)-0.5);
    end
    
    J2 = zeros(size(J));
    JID2 = sub2ind(size(J2),y,x);
    
    drMin = min(dr,[],2);
    drMin = drMin*ones(1,nFound);
    [L0,~] = find((dr==drMin)');
    
    %Ensure that each pixel is assigned to only one of the new objects by
    %adding noise.
    while numel(L0)~=numel(JID2)
        dr = dr + 1e-10*(rand(size(dx,1),nFound)-0.5);
        drMin = min(dr,[],2);
        drMin = drMin*ones(1,nFound);
        [L0,~] = find((dr==drMin)');
    end
    
    J2(JID2) = L0;
        
    %Get the boundaries of the newly separated objects
    B = cell(max(L0),1);
    for j = 1:max(L0)
        Bj = bwboundaries(J2==j,'noholes');
        if isempty(Bj)
            %If all points
            tf = 0;
            return
        else
            B{j} = Bj{1}(:,[2 1]);
        end
    end
    
    cellArea = zeros(1,numel(B));
    for j = 1:numel(B)
        A = getCellPropsMP((B{j}-1)/resizeScale);
        cellArea(j) = A.area;
    end
    if any(cellArea<cellAreaMin)
        tf = 0;
        return
    end
    
    %UPDATE CELLS1
    %Clear all fields referring to the old boundary
    chLabel = CELLS1.channel(cellID);
    
    CELLS1.boundary(cellID) = [];
    CELLS1.area(cellID) = [];
    CELLS1.centroid(cellID,:) = [];
    CELLS1.channel(cellID) = [];
    CELLS1.rank(cellID) = []; %**nb. rank will be updated elsewhere (doLabel loop)
    
    if cellID<=numel(CELLS1.label)
        CELLS1.label(cellID) = [];
    end
    
    if cellID<=numel(CELLS1.lineID);
        CELLS1.lineID(cellID) = [];
    end
    
    nCellsLeft = numel(CELLS1.area);
    
    %Input new properties
    for nf = 1:numel(B)
        CELLS1.boundary{nCellsLeft+1} = (B{nf}-1)/resizeScale + ones(size(B{nf},1),1)*bb([3 1]);
        CELLS1.area(nCellsLeft+1) = cellArea(nf);
        CELLS1.centroid(nCellsLeft+1,:) = mean(CELLS1.boundary{nCellsLeft+1});
        CELLS1.channel(nCellsLeft+1) = chLabel;
        nCellsLeft = nCellsLeft + 1;
    end
        
    tf = 1; %Resegmentation occurred and changed the number of cells found.
end

end %function resegmentMP



%%
function C = rmBorderCellsMP(C)
%C = RMBORDERCELLSMP(C).
%Filter to remove cells that are touching the bounding boxes.

deleteMe = false(numel(C.CELLS.boundary),1);
for rr = 1:numel(C.CELLS.boundary)

    xx = C.CELLS.boundary{rr}(:,1);
    yy = C.CELLS.boundary{rr}(:,2);
    
    bb = C.CHANNELS.boundingBox(C.CELLS.channel(rr),:);
    
    tf1 = any(xx<=(bb(3)+1));
    tf2 = any(xx>=(bb(4)-1));
    tf3 = any(yy<=(bb(1)+1));
    tf4 = any(yy>=(bb(2)-1));
    
    if tf1 || tf2 || tf3 || tf4
        deleteMe(rr) = true;
    end
end

%Remove corresponding cells
C.CELLS.boundary(deleteMe) = [];
C.CELLS.area(deleteMe) = [];
C.CELLS.centroid(deleteMe,:) = [];
C.CELLS.channel(deleteMe) = [];

%Check if any channels should be deleted altogether
channelsBefore = 1:numel(C.CHANNELS.threshold);
channelsLeft = unique(C.CELLS.channel);

if not(isequal(channelsBefore,channelsLeft))
    deleteMe = ~ismember(channelsBefore,channelsLeft);
    C.CHANNELS.boundingBox(deleteMe,:) = [];
    C.CHANNELS.threshold(deleteMe) = [];
    C.CHANNELS.rawIntensityRange(deleteMe,:) = [];
    C.CHANNELS.QC(deleteMe,:) = [];
    channel0 = C.CELLS.channel;

    deleteMe = find(deleteMe);
    deleteMe = reshape(deleteMe,1,numel(deleteMe));
    for rr = deleteMe
        cIDk = channel0>rr;
        C.CELLS.channel(cIDk) = C.CELLS.channel(cIDk) - 1;    
    end
end


end %function rmBorderCellsMP




end %MAIN BODY