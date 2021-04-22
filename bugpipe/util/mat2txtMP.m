function stats = mat2txtMP(D)
%stats = MAT2TXTMP(D). Writes line data to TXT. (2016-10-31 data format). D should
%be of the getDirsMP form.

lineLS = fullfile(D.lineDir,'lineDirLS.mat');
if ~exist(lineLS,'file')
    L = lineDirListMP(D,1,1);
else
    L = load(lineLS);
    L = L.L;
end

clearDirsMP(D,{'textDir','textLineDir','flatDir'})

stats = {};
stats.Ntot = 0;
stats.Nwhole = 0;
stats.Nroots = 0;
stats.Nleaves = 0;
stats.frameDrops = 0;

xDeadEnd = xDeadEndMP(D);

tic;

headersFormat = headersFormatMP;
headersFormatFlat = headersFormatFlatMP;
B = zeros(1e6,numel(headersFlat));
bc = 1;

for f = 1:numel(L)
    for k = 1:numel(L{f}.ls)
        lineData = load(L{f}.ls{k});
        
        %Check that at least one whole cell cycle was observed for this
        %cell line
        isParent = lineData.cellLineFamilyMap(lineData.cellLineFamilyMap~=0);
        isParent = unique(isParent);
        isDaughter = find(lineData.cellLineFamilyMap~=0);
        isBoth = isParent(ismember(isParent,isDaughter));
        if sum(isBoth)==0
            continue
        else
            %Write whole cell cycles, and their direct parents and daughters
            nCellsMax = numel(lineData.cellLineFamilyMap);
            cellWritten = false(1,nCellsMax);
            A = zeros(nCellsMax*20,numel(headersLine));
            ac = 1;
            for j = isBoth
                if cellWritten(j)
                    continue
                end
            
                %LOAD CELL
                cellData = load(fullfile(D.lineDir,sprintf('fov%.2d',L{f}.fov),sprintf('line%s',lineData.cellLine{j}(4:7)),strcat(lineData.cellLine{j},'.mat')));
                cellData = cellData.cellData;
                
                nPts = numel(cellData.majorAxis);
                nFrames = diff(cellData.frames)+1;
                
                %Frame was dropped
                if nPts~=nFrames
                    stats.frameDrops = stats.frameDrops + nFrames-nPts;
                end
                
                %Store LINE data
                A0 = fillAMP(cellData);
                A(ac:ac+size(A0,1)-1,:) = A0;
                
                stats.Ntot = stats.Ntot + 1;
                stats.Nwhole = stats.Nwhole + 1;
                cellWritten(j) = true;
                ac = ac+size(A0,1);
                
                %Store FLAT data
                B(bc,:) = fillBMP(cellData);
                bc = bc + 1;
                
                %CHECK PARENT
                if any(isBoth==lineData.cellLineFamilyMap(j))
                    %Parent will be written through j loop
                elseif cellWritten(lineData.cellLineFamilyMap(j))
                    %Parent was already written
                else
                    %Load parent data
                    cellData = load(fullfile(D.lineDir,sprintf('fov%.2d',L{f}.fov),...
                        sprintf('line%s',lineData.cellLine{j}(4:7)),strcat(lineData.cellLine{lineData.cellLineFamilyMap(j)},'.mat')));
                    cellData = cellData.cellData;

                    nPts = numel(cellData.majorAxis);
                    nFrames = diff(cellData.frames)+1;

                    %Frame was dropped
                    if nPts~=nFrames
                        stats.frameDrops = stats.frameDrops + nFrames-nPts;
                    end
                    
                    %Store LINE data
                    A0 = fillAMP(cellData);
                    A(ac:ac+size(A0,1)-1,:) = A0;
                                    
                    cellWritten(lineData.cellLineFamilyMap(j)) = true;
                    stats.Ntot = stats.Ntot + 1;
                    stats.Nroots = stats.Nroots + 1;
                    ac = ac+size(A0,1);
                    
                    %Store FLAT data
                    B(bc,:) = fillBMP(cellData);
                    bc = bc + 1;
                end
                
                %CHECK DAUGHTERS
                daughters = find(lineData.cellLineFamilyMap==j);
                for d = 1:numel(daughters)
                    if any(isBoth==daughters(d))
                        %Daughter will be written through j loop
                    elseif cellWritten(daughters(d))
                        %Daughter was already written
                    else
                        %Load daughter data
                        cellData = load(fullfile(D.lineDir,sprintf('fov%.2d',L{f}.fov),...
                            sprintf('line%s',lineData.cellLine{j}(4:7)),strcat(lineData.cellLine{daughters(d)},'.mat')));
                        cellData = cellData.cellData;
                        
                        nPts = numel(cellData.majorAxis);
                        nFrames = diff(cellData.frames)+1;
                        
                        %Frame was dropped
                        if nPts~=nFrames
                            stats.frameDrops = stats.frameDrops + nFrames-nPts;
                        end
                        
                        A0 = fillAMP(cellData);
                        A(ac:ac+size(A0,1)-1,:) = A0;
                
                        cellWritten(lineData.cellLineFamilyMap(j)) = true;
                        stats.Ntot = stats.Ntot + 1;
                        stats.Nleaves = stats.Nleaves + 1;
                        ac = ac+size(A0,1);
                        
                        %Store FLAT data
                        B(bc,:) = fillBMP(cellData);
                        bc = bc + 1;
                    end
                end                
            end
            
            A = A(1:ac-1,:);
                        
            %Write to TXT
            lineTXT = fullfile(D.textDir,'lineages',[lineData.cellLine{1}(1:7) '.txt']);
            fidTXT = fopen(lineTXT,'w');
            fprintf(fidTXT,headersFormat,A');
            fclose(fidTXT);
            
            %Save to MAT
            save(fullfile(D.flatDir,[lineData.cellLine{1}(1:7) '.mat']),'A')
        end        
    end
end
B = B(1:bc-1,:);
%Write to TXT
flatTXT = fullfile(D.textDir,[D.baseDir '_flat.txt']);
fidTXT = fopen(flatTXT,'w');
fprintf(fidTXT,headersFormatFlat,B');
fclose(fidTXT);

%Save to MAT
save(fullfile(D.flatDir,[D.baseDir '_flat.mat']),'B','stats','headersLine','headersFlat');

function headersFormat = headersFormatMP
    headersLine = {'fov','line','cellID','parentID','frame','t','L','w','area','Itot','orientation','xCM','yCM','xCMdead','xCMweight','yCMweight','QC'};
    headersFormat0 = cell(numel(headersLine),1);
    headersFormat0{indMP1('fov')} = '%.2d';
    headersFormat0{indMP1('line')} = '%.4d';
    headersFormat0{indMP1('cellID')} = '%d';
    headersFormat0{indMP1('parentID')} = '%d';
    headersFormat0{indMP1('frame')} = '%d';
    theRest = headersLine(~ismember(headersLine,{'fov','line','frame'}));
    for h = 1:numel(theRest)
        headersFormat0{indMP1(theRest{h})} = '%7.7f';
    end
    
    headersFormat = headersFormat0{1};
    for h = 2:numel(headersLine)
         headersFormat = [headersFormat '\t' headersFormat0{h}];
    end
    
    headersFormat = [headersFormat '\n'];
end

function headersFormatFlat = headersFormatFlatMP
    headersFlat = {'fov','line','cellID','parentID','frame0','t0','divT','alphaL','alphaV','L0','L1','w0','w1','area0','area1',...
        'V0','V1','L0fit','L1fit','V0fit','V1fit','Itot0','Itot1','orientation0','orientation1','xCM0','yCM0','xCM1','yCM1',...
        'xCMdead0','xCMdead1','xCMweight0','yCMweight0','xCMweight1','yCMweight1','QC'};
    headersFormat0 = cell(numel(headersFlat),1);
    headersFormat0{indMP2('fov')} = '%.2d';
    headersFormat0{indMP2('line')} = '%.4d';
    headersFormat0{indMP2('cellID')} = '%d';
    headersFormat0{indMP2('parentID')} = '%d';
    headersFormat0{indMP2('frame0')} = '%d';
    theRest = headersFlat(~ismember(headersFlat,{'fov','line','frame'}));
    for h = 1:numel(theRest)
        headersFormat0{indMP2(theRest{h})} = '%7.7f';
    end
    
    headersFormatFlat = headersFormat0{1};
    for h = 2:numel(headersFlat)
         headersFormatFlat = [headersFormatFlat '\t' headersFormat0{h}];
    end
    
    headersFormatFlat = [headersFormatFlat '\n'];
end

function ID = indMP1(prop)
    ID = strcmp(prop,headersLine);
end

function ID = indMP2(prop)
    ID = strcmp(prop,headersFlat);
end

function A0 = fillAMP(cellData)
    fov = str2double(cellData.label(1:2));
    pickMe = true(size(cellData.majorAxis));
    xDead = xDeadEnd.x(xDeadEnd.fov==fov & xDeadEnd.t >= cellData.frames(1) & xDeadEnd.t <= cellData.frames(2));
    if numel(xDead)~=numel(cellData.majorAxis)
        xDead0 = xDead(xDead~=0); %If a frame was dropped
        if numel(xDead0)==numel(cellData.majorAxis)
            xDead = xDead0;
        elseif numel(cellData.majorAxis)>numel(xDead);
            %2016-11-08: in all datasets so far (33 in dataBook), this error
            %has only occurred ONCE. apparently due to double-indexing, i.e.
            %two frames were each registered twice so that there were two
            %extra data points in each property BUT the frame span remained
            %the same. This is obvious because time acquired repeats.
            %SEE ..\Data_err\doubleFrames
            dt = diff(cellData.tAcq);
            pickMe = [true(1); logical(dt)];
        end
    end
    
    A0 = zeros(sum(pickMe),numel(headersLine));
    
    pix2um = 0.1067;
    A0(:,indMP1('fov')) = fov;
    A0(:,indMP1('line')) = str2double(cellData.label(4:7));
    A0(:,indMP1('cellID')) = cellData.lineID;
    A0(:,indMP1('parentID')) = cellData.parent;
    A0(:,indMP1('frame')) = cellData.frames(1):cellData.frames(1)+sum(pickMe)-1;
    A0(:,indMP1('t')) = cellData.tAcq(pickMe);
    A0(:,indMP1('L')) = cellData.majorAxis(pickMe)*pix2um;
    A0(:,indMP1('w')) = cellData.minorAxis(pickMe)*pix2um;
    A0(:,indMP1('area')) = cellData.area(pickMe)*(pix2um^2);
    A0(:,indMP1('Itot')) = cellData.intensityMeanSD(pickMe);
    A0(:,indMP1('orientation')) = cellData.orientation(pickMe);
    A0(:,indMP1('xCM')) = cellData.regionpropsCentroid(pickMe,1)*pix2um;
    A0(:,indMP1('yCM')) = cellData.regionpropsCentroid(pickMe,2)*pix2um;
    A0(:,indMP1('xCMweight')) = cellData.weightedCentroid(pickMe,1)*pix2um;
    A0(:,indMP1('yCMweight')) = cellData.weightedCentroid(pickMe,2)*pix2um;
    A0(:,indMP1('xCMdead')) = abs(cellData.regionpropsCentroid(pickMe,1) - xDead)*pix2um;
    A0(:,indMP1('QC')) = cellData.QC(pickMe);
end

function B0 = fillBMP(cellData)
    B0 = zeros(1,numel(headersFlat));
    fov = str2double(cellData.label(1:2));
    
    vol = @(L,w) 4/3*pi*(w/2).^3 + pi*(w/2).^2.*(L-w);
    xDead = xDeadEnd.x(xDeadEnd.fov==fov & xDeadEnd.t >= cellData.frames(1) & xDeadEnd.t <= cellData.frames(2));
    xDead = xDead(xDead~=0);
    
    pix2um = 0.1067;
    B0(indMP2('fov')) = fov;
    B0(indMP2('line')) = str2double(cellData.label(4:7));
    B0(indMP2('cellID')) = cellData.lineID;
    B0(indMP2('parentID')) = cellData.parent;
    B0(indMP2('frame0')) = cellData.frames(1);
    B0(indMP2('t0')) = cellData.tAcq(1);
    B0(indMP2('divT')) = cellData.tAcq(end) - cellData.tAcq(1) + xDeadEnd.dt;
    B0(indMP2('L0')) = cellData.majorAxis(1)*pix2um;
    B0(indMP2('L1')) = cellData.majorAxis(end)*pix2um;
    B0(indMP2('w0')) = cellData.minorAxis(1)*pix2um;
    B0(indMP2('w1')) = cellData.minorAxis(end)*pix2um;
    B0(indMP2('area0')) = cellData.area(1)*(pix2um^2);
    B0(indMP2('area1')) = cellData.area(end)*(pix2um^2);
    B0(indMP2('V0')) = vol(cellData.majorAxis(1),cellData.minorAxis(1))*(pix2um^3);
    B0(indMP2('V1')) = vol(cellData.majorAxis(end),cellData.minorAxis(end))*(pix2um^3);
    B0(indMP2('Itot0')) = cellData.intensityMeanSD(1);
    B0(indMP2('Itot1')) = cellData.intensityMeanSD(end);
    B0(indMP2('orientation0')) = cellData.orientation(1);
    B0(indMP2('orientation1')) = cellData.orientation(end);
    B0(indMP2('xCM0')) = cellData.regionpropsCentroid(1,1)*pix2um;
    B0(indMP2('yCM0')) = cellData.regionpropsCentroid(1,2)*pix2um;
    B0(indMP2('xCM1')) = cellData.regionpropsCentroid(end,1)*pix2um;
    B0(indMP2('yCM1')) = cellData.regionpropsCentroid(end,2)*pix2um;
    B0(indMP2('xCMweight0')) = cellData.weightedCentroid(1,1)*pix2um;
    B0(indMP2('yCMweight0')) = cellData.weightedCentroid(1,2)*pix2um;
    B0(indMP2('xCMdead0')) = abs(cellData.regionpropsCentroid(1,1) - xDead(1))*pix2um;
    B0(indMP2('xCMdead1')) = abs(cellData.regionpropsCentroid(end,1) - xDead(end))*pix2um;
    B0(indMP2('xCMweight1')) = cellData.weightedCentroid(end,1)*pix2um;
    B0(indMP2('yCMweight1')) = cellData.weightedCentroid(end,2)*pix2um;
    B0(indMP2('QC')) = mean(cellData.QC);
    
    %Get growth rates and fitted birth/division values
    v = vol(cellData.majorAxis,cellData.minorAxis)*(pix2um^3);
    l = cellData.majorAxis*pix2um;
    t = cellData.tAcq - cellData.tAcq(1);
    
    if numel(t)<2
        B0(indMP2('alphaV')) = NaN;
        B0(indMP2('V0fit')) = v;
        B0(indMP2('V1fit')) = v;
        
        B0(indMP2('alphaL')) = NaN;
        B0(indMP2('L0fit')) = l;
        B0(indMP2('L1fit')) = l;
    else
        f0 = fit(t(:),v(:),'exp1','startpoint',[1 1/60*log(2)]);
        B0(indMP2('alphaV')) = f0.b;
        B0(indMP2('V0fit')) = f0.a*exp(f0.b*(-xDeadEnd.dt/2));
        B0(indMP2('V1fit')) = f0.a*exp(f0.b*(cellData.tAcq(end)+xDeadEnd.dt/2));

        f0 = fit(t(:),l(:),'exp1','startpoint',[1 1/60*log(2)]);
        B0(indMP2('alphaL')) = f0.b;
        B0(indMP2('L0fit')) = f0.a*exp(f0.b*(-xDeadEnd.dt/2));
        B0(indMP2('L1fit')) = f0.a*exp(f0.b*(cellData.tAcq(end)+xDeadEnd.dt/2));
    end
end

end