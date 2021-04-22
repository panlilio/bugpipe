function P = segmentationParamsMP(leftRight)
%P = SEGMENTATIONPARAMSMP.
%Initializes segmentation parameters to be passed to segmentationMP.

P = {};
P.resizeScale = 3;      %Scale for imresize to find smooth boundaries.
P.se = strel('disk',3); %Structuring element for imdilate
P.cellAreaMin = 30;

if nargin==0
    P.deadEnd = 'left';
else
    P.deadEnd = leftRight;
end

end
