function [x0,y0] = constructBoundaryMP(xx,yy)
%[x,y] = constructBoundaryMP(xx,yy).
%Outputs boundary points in whole pixels so that imfill can be performed.

%Catch erroneous objects and do not try to reconstruct a boundary.
%Note: as of 2016-10-31 this error has only occurred in unusable frames
%that would be neglected during analysis. e.g. See end frames of 
%2015-02-23 dataset, where a real segmentation cannot be obtained.
if (max(xx)-min(xx))>512 || (max(yy)-min(yy))>512 
    x0 = round(xx(:));
    y0 = round(yy(:));
    return
end

%NORMAL PROCEDURE: bridge pixel gaps between points. 
x0 = xx(:);
y0 = yy(:);

x0 = round(x0);
y0 = round(y0);

x0 = [x0(:); x0(1)];
y0 = [y0(:); y0(1)];

x = x0;
y = y0;

dr = diff(x0).^2 + diff(y0).^2;

tooFar = find(dr>2);
while any(tooFar)
    tID = 0*dr;
    tID(tooFar+1) = 1;
    tID = cumsum(tID);
    tID = tID + 1:numel(dr);
    
    tooFarMean = [tooFar'; tooFar'+1];
    x = x0(tID);
    y = y0(tID);
    
    xMean = mean(x0(tooFarMean));
    yMean = mean(y0(tooFarMean));
    
    x(tooFar) = xMean;
    y(tooFar) = yMean;
    
    x = round(x);
    y = round(y);
    
    x0 = [x; x(1)];
    y0 = [y; y(1)];
    dr = diff(x0).^2 + diff(y0).^2;
    
    tooFar = find(dr>2);
end

% Eliminate any redundant points
xy = [x(:) y(:)];
% xy = unique(xy,'rows');
x0 = [xy(:,1); xy(1,1)];
y0 = [xy(:,2); xy(1,2)];