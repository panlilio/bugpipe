function A = getCellPropsMP(boundaryPts,I,tAcq)
%A = GETCELLPROPSMP(boundaryPts,I,tAcq). Calculates properties of the cell
%with boundary points boundaryPts. Determines intensity values using image
%I and records the time that the cell was observed, tAcq.

A = {};

%Calculate area
xx = boundaryPts(:,1);
yy = boundaryPts(:,2);

[xx,yy] = constructBoundaryMP(xx,yy);

xx = min(512,max(1,xx));
yy = min(512,max(1,yy)); 

xx0 = min(xx);
yy0 = min(yy);

xxSmall = xx - xx0 + 1;
yySmall = yy - yy0 + 1;

cmyx = round(mean([yySmall xxSmall],1));
J = false(max(yySmall),max(xxSmall));
JID = sub2ind(size(J),yySmall,xxSmall);
J(JID) = true;

J = imfill(J,cmyx,4);

%COMMENTED OUT. The following method should be equivalent and faster than
%imfill but it seems to produce different results and cell tracking
%changes (increased number of very short cell divisions).
% L = bwlabel(not(J)); 
% for n = 1:max(L(:))
%     [r,c] = find(L==n);
%     if any(r==size(J,1)) || any(c==size(J,2))
%         J(L==n) = false;
%     else
%         J(L==n) = true;
%     end
% end

A.area = sum(J(:));

if nargin==1 || isempty(I)
    return
end

%Calculate rest of properties
Ismall = I(max(1,yy0):max(1,yy0)+size(J,1)-1,max(1,xx0):max(1,xx0)+size(J,2)-1); %cropped image
Ismall(not(J)) = 0; 

%Intensity
intensityVals = Ismall(logical(J));
meanI = sum(intensityVals);

A.intensityRange = [min(intensityVals) max(intensityVals)];
A.intensityMeanSD = meanI;

%PCA to get the ellipse axes and orientation.
[coeff,score] = pca(boundaryPts);
p1 = score(:,1);
p2 = score(:,2);

dp1 = p1*ones(1,numel(p1));
dp1 = dp1 - dp1';
dp1 = abs(dp1);
dp2 = p2*ones(1,numel(p2));
dp2 = dp2 - dp2';
dp2 = abs(dp2);

A.majorAxis = max(dp1(:));
A.minorAxis = max(dp2(:));
A.orientation = atan(coeff(2)/coeff(1));

[row,col,val] = find(Ismall);
A.weightedCentroid = sum([col.*val row.*val])./sum(val) + [xx0 yy0] - 1;
A.regionpropsCentroid = mean([col row]) + [xx0 yy0] - 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%20170401 ADDED: surface area, volume, and alternative width calculations
%Translate to centroid
xx = boundaryPts(:,1) - mean(boundaryPts(:,1));
yy = boundaryPts(:,2) - mean(boundaryPts(:,2));

%Rotate to principal axes
R = [cos(A.orientation), sin(A.orientation); -sin(A.orientation) cos(A.orientation)];
xy = [reshape(xx,1,numel(xx)); reshape(yy,1,numel(yy))];
xy = R*xy;
xy = xy';
xy = [xy; xy(1,:)];

%Cycle indices so that positive and negative y segments are continuous
pos = xy(:,2)>=0;
pos = [pos(1); diff(pos)];
f = find(pos==1,1,'first');
xy = circshift(xy,1-f,1);

%Calculate surface area as the mean of rotations of positive and negative
%perimeters
pos = xy(:,2)>=0;
neg = xy(:,2)<=0;

dxy = [diff(xy); xy(end,:)-xy(1,:)];
dr = sqrt(dxy(:,1).^2 + dxy(:,2).^2);

SApos = abs(2*pi*sum(xy(pos,2).*dr(pos)));
SAneg = abs(2*pi*sum(xy(neg,2).*dr(neg)));

A.surfaceArea = mean([SApos SAneg]);

%Calculate volume as the mean of rotations of positive and negative
%perimeters
dx = dxy(:,1);

Vpos = abs(pi*sum(xy(pos,2).^2.*dx(pos)));
Vneg = abs(pi*sum(xy(neg,2).^2.*dx(neg)));

A.volume =  mean([Vpos Vneg]);

%ALTERNATIVE WIDTH CALCULATIONS
%As the V/SA*4 (according to cylinder scaling)
A.widthCyl = A.volume./A.surfaceArea*4;
A.lengthCyl = (A.surfaceArea.^2)./A.volume./(4*pi)+A.widthCyl;

%Assuming the area is a projected spherocyclinder
A.widthSphcyl = 2*(A.majorAxis+1 - sqrt((A.majorAxis+1).^2-(4-pi).*A.area))./(4-pi) - 1;

%As the mean width, of the middle of the cell (by angle)
xypos = xy(pos,:);
xyneg = xy(neg,:);

try
    %Defining the "middle" of the cell sometimes fails due to bad
    %boundaries, i.e. theta is badly behaved
    
    theta = abs(atan(xypos(:,2)./xypos(:,1))); %angular coord of positive cell boundary
    thetaRange = [find(theta>pi/4,1,'first') find(theta>pi/4,1,'last')]; %points between pi/4 and 3pi/3. Note this only works because x is ordered.
    xypos = xypos(thetaRange(1):thetaRange(2),:); %limit the positive perimeter to these points
    
    %Find the distance across the cell by minimizing the x-distance
    deltaX = ones(size(xyneg,1),1)*(xypos(:,1))' - xyneg(:,1)*ones(1,size(xypos,1));
    deltaR = deltaX.^2;

    [~,minID] = min(deltaR,[],1);
    xypair1 = xypos;
    xypair2 = xyneg(minID,:);

    deltaR = xypair1 - xypair2;
    deltaR = sqrt(deltaR(:,1).^2 + deltaR(:,2).^2);

    A.widthMean = mean(deltaR);

catch
    A.widthMean = NaN;
end

%20170401 END.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Determine the quality of segmentation: check cross sections for
%double-loaded channels
dI = diff(Ismall);
dI = sign(dI);
dI = diff(dI);
[row,col] = find(dI==-2);

nLocalMax = 0;
uncol = unique(col);
for j = reshape(uncol,1,numel(uncol))
    r = find(J(:,j));
    fMax = row(col==j);
    fMax = fMax(fMax~=min(r) & fMax~=max(r));
    if any(diff(fMax)>2)
        nLocalMax = nLocalMax + 1;
    end
end
A.QC = nLocalMax/numel(uncol);

if nargin==2
    A.tAcq = NaN;
    return
end

%Record time that the cell was observed
A.tAcq = tAcq;