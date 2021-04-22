function [CELLS,CHANNELS] = segmentationMP(I,params)
%[CELLS,CHANNELS] = segmentationMP(I). Cell segmentation of the intensity matrix I
%from mother machine data, output CELLS is a structured array

%Segmentation parameters
resizeScale = params.resizeScale;
se = params.se; 
cellAreaMin = params.cellAreaMin;

%Suppress common warning
warning('off','MATLAB:rankDeficientMatrix')

%Initialize variables
CELLS = {};
CELLS.boundary = {};
CELLS.area = 0;
CELLS.centroid = [0 0];
CELLS.channel = 0;
cellCount = 1;

CHANNELS = {};

%Get channels
[bb,QC] = getChannelsMP(I);
nChannels = size(bb,1);
if nChannels > 11 || nChannels==0
    fprintf('\t ERROR: No channels or too many found--corrupt image? \n')
    return
end

CHANNELS.boundingBox = bb;
CHANNELS.QC = QC;
CHANNELS.threshold = zeros(1,nChannels);
CHANNELS.rawIntensityRange = zeros(nChannels,2);

%Get cells in each channel
for k = 1:nChannels
    Ik1 = I(bb(k,1):bb(k,2),bb(k,3):bb(k,4));
    Ik1 = medfilt2(Ik1);
    Ik2 = imresize(Ik1,resizeScale,'bilinear');
    Ik2 = mat2gray(Ik2);
    
    %Threshold image to binary
    tOtsu = graythresh(Ik2);
    lowI = Ik2(Ik2(:)<tOtsu);
    sig = sqrt(1/numel(lowI)*sum(lowI.^2)-(1/numel(lowI)*sum(lowI)).^2);
    T = tOtsu + 1*sig;
    BW = Ik2>=T;
    BW = imdilate(BW,se);
    
    %Record channel image properties
    CHANNELS.rawIntensityRange(k,:) = [min(Ik1(:)) max(Ik1(:))];
    CHANNELS.threshold(k) = T;
    
    %Get boundaries of binary image
    B = bwboundaries(BW,'noholes');
    for j = 1:numel(B)
        CELL0 = getCellsMP(B{j});
        CELL0 = CELL0(not(cellfun(@isempty,CELL0)));
        
        %Filter by area
        for nc = 1:numel(CELL0)
            A = getCellPropsMP((CELL0{nc}-1)/resizeScale);
            cellArea = A.area;
            if cellArea >= cellAreaMin
                CELLS.boundary{cellCount} = (CELL0{nc}-1)/resizeScale + ones(size(CELL0{nc},1),1)*bb(k,[3 1]);
                CELLS.area(cellCount) = cellArea;
                CELLS.centroid(cellCount,:) = mean(CELLS.boundary{cellCount});
                CELLS.channel(cellCount) = k;
                cellCount = cellCount + 1;
            end
        end
    end %processing object in each boundary    
end %channels


clearvars -except CELLS CHANNELS


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ADDITIONAL FUNCTIONS CALLED BY MAIN BODY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [bb,QC] = getChannelsMP(I)
%[bb,QC] = GETCHANNELSMP(I). Finds channels based on high intentsity
%regions of the image I. Outputs the pixel indices of each channel.
%bb is an n by 4 matrix  such that 
%bb(k,:) = [rowMin rowMax colMin colMax] of the kth channel. Also outputs
%a quality control matrix QC which contains the number of local maxima in the
%channel cross sections to detect either empty channels (very noisy) or
%double-loaded channels (will contain multiple maxima per cross-section).

%Resample image at lower spatial frequency to determine approximate channel edges
sampleRate = 5;
Ismall = I(1:sampleRate:end,1:sampleRate:end);

Ix = max(Ismall,[],1); Ix = Ix(:);
Iy1 = mean(Ismall,2); Iy1 = Iy1(:);
Iy2 = smooth(Iy1,4*sampleRate+1);
Iy = Iy1 - Iy2;
Iy = max(0,Iy);
Iy = logical(Iy);

%Find the bounding box of each channel. %Consider changing these values.
epsX = 10; %pixel cushion in x-direction to expand found bounding box
epsY = 10; %-------"-------- y-direction ------------"---------------

%x-bound
IxThreshold = graythresh(mat2gray(Ix)).*(max(Ix)-min(Ix))+min(Ix);
Ix = Ix>=IxThreshold;
dIx = [sign(diff(Ix)); 0];
x1 = find(dIx>0,1,'first'); if isempty(x1), x1 = 1; end
x2 = find(dIx<0,1,'last'); if isempty(x2), x2 = floor(512/sampleRate); end

dx = ([x1 x2]-1)*sampleRate+1; %Translate back to full image index
dx(1) = max(1,dx(1)-epsX);
dx(2) = min(512,dx(2)+epsX);

%y-bound
dIy = [sign(diff(Iy)); 0];
y1 = find(dIy>0); if isempty(y1), y1 = 1; end;
y2 = find(dIy<0); if isempty(y2), y2 = floor(512/sampleRate); end

if y2(1) < y1(1); y1 = [1; y1]; end

m = min(numel(y1),numel(y2));
y1 = y1(1:m);
y2 = y2(1:m);

dy = ([y1 y2]-1)*sampleRate+1; 
dy(:,1) = max(1,dy(:,1)-epsY);
dy(:,2) = min(512,dy(:,2)+epsY);

dyNew = dy*0;
dyNew(1,:) = dy(1,:);
c = 1;
for k1 = 2:m
    if dyNew(c,2)>dy(k1,1)
        dyNew(c,2) = dy(k1,2);
    else
        c = c+1;
        dyNew(c,:) = dy(k1,:);
    end
end
dy = dyNew(1:c,:);

%Check that there are cells present in each interval of the image, store
%the bounding box and cross-section quality control metric if so.
bb = 0*[dyNew dyNew];
QC = zeros(size(bb,1),2);
c = 1;
for k1 = 1:size(dy,1)
    Ik = I(dy(k1,1):dy(k1,2),dx(1):dx(2));
    if sum(Ik(:)>500)<=10 %Need to modify this for gyrA
        %No high intensity objects: discard channel.
    else
        %Check y-cross section for well-defined objects (i.e. approximately
        %gaussian intensity distributions with a single peak), save results
        %in quality control matrix QC
        Ik = mat2gray(Ik);
        thresh = graythresh(Ik);
        Ik(Ik<=thresh) = 0;
        dIk = diff(Ik,1,1);
        dIk = sign(dIk);
        ddIk = diff(dIk);
        [row,col] = find(ddIk==-2);
        validC = unique(col);
        validC = reshape(validC,1,numel(validC));
        qc = [0 0];
        for v = validC
            if any(diff(row(col==v))>2)
                %There are multiple maxima spaced more than two pixels apart
                qc(1) = qc(1) + 1;
            end
        end
        %Total number of cross-sections with any values above threshold
        qc(2) = numel(validC);
                
%         if qc(1)>80
            %Too many multiple maxima: likely an empty channel.
%         else
            %Record boundary box points and quality control parameters
            QC(c,:) = qc;
            bb(c,:) = [dy(k1,1) dy(k1,2) dx(1) dx(2)];
            c = c + 1;
%         end
    end
end
bb = bb(1:c-1,:);
QC = QC(1:c-1,:);

end %GETCHANNELSMP function




end