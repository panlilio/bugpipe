function I = getCellIntensityMP(boundaryPts,I)
%I = GETCELLINTENSITYMP. Retrieves the intensity values of the cell
%with boundary points boundaryPts. Code for reconstruction of cell is taken
%from getCellPropsMP--therefore update this as necessary.

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

%Calculate rest of properties
Ismall = I(max(1,yy0):max(1,yy0)+size(J,1)-1,max(1,xx0):max(1,xx0)+size(J,2)-1); %cropped image
Ismall(not(J)) = 0; 

%Intensity
I = Ismall(logical(J));
