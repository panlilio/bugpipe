classdef dataMP2 < handle
    properties 
        dirList
        strain
        orientation
        switchTime
        uniformLight
        frameLim
        media
        comments
        flatFile
        lineFiles
        headersLine
        headersFlat
        a
        b
        mat2txtstats
        filterIDFlat
        cellFilterFlat
        filterIDLine
        cellFilterLine
        smoothingMethod
        flatData
        imList
        parentList
        uniqueID
        aIDmap
        diffaID
        currShift
        storedFilters
        figSaveDir
    end
    
    methods
        
    function obj = dataMP2(D)
    %OBJ = DATAMP2(D)
    %Initialize dataMP2 object with data directory structure D, output from
    %getDirsMP or the string of format 'baseDir'
    if ~isstruct(D)
        D = getDirsMP_SSD(D);
    end
    
        dataInfo = lookupDataMP(D.baseDir);
        
        obj.dirList = D;
        obj.strain = dataInfo.STRAIN;
        obj.orientation = dataInfo.ORIENTATION;
        obj.switchTime = dataInfo.SWITCH;
        obj.uniformLight = dataInfo.UNIFORMLIGHT;
        obj.frameLim = dataInfo.FRAMELIM;
        obj.media = dataInfo.MEDIA;
        obj.comments = dataInfo.COMMENTS;
        
        flatData = load(fullfile(D.flatDir,[D.baseDir '_flat.mat']));
        obj.flatData = flatData;
        
        obj.flatFile = fullfile(D.flatDir,[D.baseDir '_flat.mat']);
        obj.headersLine = obj.flatData.headersLine;
        obj.headersFlat = obj.flatData.headersFlat;
        obj.mat2txtstats = obj.flatData.stats;
        
        obj.a = containers.Map(obj.headersLine,1:numel(obj.headersLine));
        obj.b = containers.Map(obj.headersFlat,1:numel(obj.headersFlat));
                
        obj = obj.initFilter;
        
        %Smoothing method to use for e.g. area, length, width, volume
        obj.smoothingMethod = struct('dosmooth',false,'A',struct,'B',struct);
        
        if ispc
            L = ls(D.flatDir);
            L = L(3:end,:);
            L = cellstr(L);
            L = L(~cellfun(@any,strfind(L,'_flat.mat')));
        else
            L = dir(D.flatDir);
            notFlat = false(numel(L),1);
            for j = 1:numel(L)
                if ~any(strfind(L(j).name,'_flat.mat')) && numel(L(j).name)>2
                    notFlat(j) = true;
                end
            end
            L = L(notFlat);
            L = vertcat(L(:).name);
            L = cellstr(L);
        end
        obj.lineFiles = L;
                
        obj = obj.initUniqueID;
        
        obj.aIDmap = struct;
        obj.aIDmap(1).built = false;
        
        obj.diffaID = struct;
        
        obj.currShift = 0; %shift of interest, around which all times are translated, if requested 
        
        %Path to MAT with containers of stored data filters
        dataHandlingDir = which('dataMP2');
        [dataHandlingDir,~,~] = fileparts(dataHandlingDir);
        obj.storedFilters = load(fullfile(dataHandlingDir,'dataFilters.mat'));

        %Directory to save any figures
        obj.figSaveDir = '';
    end
    
    
    
    function obj = getImList(obj)
    %obj = GETIMLIST(obj)
    %Saves the output from the function imListMP: retrieves a list of all
    %png images, including blanks, and functions to look up the names of
    %png files based on field of view and time point.
    if isempty(obj.imList)
        obj.imList = imListMP(obj.dirList.dataDir);
    end
    
    end
    
    
    
    function [fov,line] = fovLineK(obj,K)
    %[fov,line] = obj.fovLineK(K)
    %Field of view and lineage number of the Kth line file in the dataset.
    
    fov = obj.lineFiles(K);
    fov = vertcat(fov{:});
    line = str2num(fov(:,4:7)); %#ok<ST2NM>
    fov = str2num(fov(:,1:2)); %#ok<ST2NM>
    
    end
    
    
    function A = loadLineK(obj,K)
    %A = loadLineK(obj,K)
    %Loads A matrix for Kth line
    
    A = load(fullfile(obj.dirList.flatDir,obj.lineFiles{K}));
    A = A.A;
    
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% FLAT: CELL FILTERS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    function obj = initFilter(obj)
    %obj = INITFILTER(obj)
    %Initialize filter for flat cell data. Default is all-pass, i.e. for
    %cell properties the allowed range is [-Inf Inf], symmetric division
    %tolerance is set to Inf, and cell cycles do not have to be whole.
    obj.filterIDFlat = true(size(obj.flatData.B,1),1);

    obj.cellFilterFlat = {};
    for k = 1:numel(obj.flatData.headersFlat)
        obj.cellFilterFlat.(obj.flatData.headersFlat{k}) = [-Inf Inf];
    end
    obj.cellFilterFlat.symmDivTolP = [-Inf Inf];
    obj.cellFilterFlat.symmDivTolC = Inf;
    obj.cellFilterFlat.wholeCC = false;
    obj.cellFilterFlat.t1 = [-Inf Inf];
    obj.cellFilterFlat.wholeCCRobust = false; %%NB: this should always be the last field
    obj.cellFilterFlat.omitShift = false; %omits cells that are born between two media (not required to be robust)
    
    %Initialize filter based on line data, note though that it is indexed
    %by cell identity in B
    obj.filterIDLine = true(size(obj.flatData.B,1),1);
    
    obj.cellFilterLine = {};
    for k = 1:numel(obj.flatData.headersLine)
        obj.cellFilterLine.(obj.flatData.headersLine{k}) = [-Inf Inf];
    end
    
    end
    
    function obj = setFlatFilter(obj,propName,propRange)
    %obj = SETFLATFILTER(obj,propName,propRange)
    %Set specific filter and update the filterIDFlat property accordingly.
    %This filter will apply when calling obj.flat_cellProps, 
    %obj.line_cellProps, and obj.line_cellPropsSingleCell.Note that 
    %propName can be a cell containing multiple properties to filter, but
    %propRange must also be a cell with corresponding indices.
    
    if nargin>1
        if ~iscell(propName), propName = {propName}; propRange = {propRange}; end
        
        %Save filter values
        for k = 1:numel(propName)
            obj.cellFilterFlat.(propName{k}) = propRange{k};
        end
    end
    
        %Recalculate logical vector
        obj.filterIDFlat = true(size(obj.flatData.B,1),1);
        filterNames = fields(obj.cellFilterFlat);
        for k = 1:numel(filterNames)
            filterRange = obj.cellFilterFlat.(filterNames{k});
            if (numel(filterRange)>1 && any(isfinite(filterRange))) || (numel(filterRange)==1 && isfinite(filterRange) && filterRange~=0)
                F = obj.getFlatFilter(filterNames{k},filterRange);
                obj.filterIDFlat = obj.filterIDFlat & F;
            end
        end
        
        obj.filterIDFlat = obj.filterIDFlat & obj.filterIDLine;
    end
    
    function F = getFlatFilter(obj,propName,propRange)
    %F = GETFLATFILTER(obj,propName,propRange)
    %F is a logical vector corresponding to indices of the flat data matrix
    %obj.flatData.B. 1 if the cell passes the specified filter, 0 if not.
    %Note that this function only handles one property filter at a time.
    
    %Check if this filter was already calculated and stored by
    %wholeCCall.m
    if isfield(obj.storedFilters,propName)
        if all(obj.storedFilters.(propName)('filterVal')==propRange) && isKey(obj.storedFilters.(propName),obj.dirList.baseDir)
            %Filter was calculated using the correct value range
            F = obj.storedFilters.(propName)(obj.dirList.baseDir);
            return
        end
    end
    
    %Otherwise perform calculation now
    if strcmp(propName,'symmDivTolP') 
        %Cell passes if its division is symmetric, using all daughters
        %that were tracked

        [parentsID,childrenID] = obj.getChildren;
        notFilament = false(size(obj.flatData.B,1),1);
        for k = 1:numel(parentsID)
            pArea = obj.flatData.B(parentsID(k),obj.b('area1'));
            cArea = obj.flatData.B(childrenID{k},obj.b('area0'));
            tf = true;

            for j = 1:numel(cArea)
                tf = tf && (cArea(j)>=propRange(1)*pArea) && (cArea(j)<=propRange(2)*pArea);
            end
            notFilament(parentsID(k)) = tf;
        end
        F = notFilament;

    elseif strcmp(propName,'symmDivTolC')
        %40/60 symmetry WRT one another
        [parentsID,childrenID] = obj.getChildren;
        notFilament = false(size(obj.flatData.B,1),1);
        for k = 1:numel(parentsID)
            cArea = obj.flatData.B(childrenID{k},obj.b('area0'));
            if numel(cArea)<2
                tf = false;
            else
                tf = min(cArea)./max(cArea);
                tf = tf>=propRange;
            end
            notFilament(parentsID(k)) = tf;
        end
        F = notFilament;

    elseif strcmp(propName,'wholeCC')
        %Only look at whole cell cycles
        F = false(size(obj.flatData.B,1),1);
        if propRange
            wholeCC = obj.getWholeCC;
            F(wholeCC) = true; 
        else
            F = true(size(obj.flatData.B,1),1);
        end

    elseif strcmp(propName,'wholeCCRobust')
        %Only looks at whole cell cycles who themselves and whose 
        %parents and children pass all existing cell filters
        F = false(size(obj.flatData.B,1),1);
        if propRange
            wholeCC = obj.getWholeCCRobust;
            F(wholeCC) = true;
        else
            F = true(size(obj.flatData.B,1),1);
        end
    elseif strcmp(propName,'omitShift')
        %Omits cells that were born before and divide after the media shift
        %that has been set by the class property currShift
        F = true(size(obj.flatData.B,1),1);
        if propRange
            t0 = obj.getPropBIndex((1:size(obj.flatData.B,1))','t0');
            t1 = obj.getPropBIndex((1:size(obj.flatData.B,1))','t1');
            F((t0<obj.currShift) & (t1>obj.currShift)) = false;
        end
    else
        try
            k = obj.getPropBIndex((1:size(obj.flatData.B,1))',propName);
            F =  (k>=propRange(1)) & (k<=propRange(2));
        catch
            warning('Did not retrieve cell property %s.\n',propName)
            F = true(size(obj.flatData.B,1),1);
        end
    end

    end
    
    
    function obj = setLineFilter(obj,propName,propRange)
    %obj = obj.setLineFilter(obj,propName,propRange)
    if ~iscell(propName), propName = {propName}; end
    if ~iscell(propRange), propRange = {propRange}; end
    
    %Save filter values
    for k = 1:numel(propName)
        obj.cellFilterLine.(propName{k}) = propRange{k};
    end
    
    %Recalculate logical vector
    obj.filterIDLine = true(size(obj.flatData.B,1),1);
    filterNames = fields(obj.cellFilterLine);
    for k = 1:numel(filterNames)
        filterRange = obj.cellFilterLine.(filterNames{k});
        if (numel(filterRange)>1 && any(isfinite(filterRange))) 
            F = obj.getLineFilter(filterNames{k},filterRange);
            obj.filterIDLine = obj.filterIDLine & F;
        end
    end
    
    %Update filterIDFlat, since this is the filter that is used to assess
    %the validity of cells during calculations and lineage building
    obj.setFlatFilter; %Updates filterIDFlat
    
    end
    
    
    function F = getLineFilter(obj,propName,propRange)
    %F = obj.getLineFilter(propName,propRange)
    %Retrieves the logical vector corresponding to indices in B of cells
    %passing the line data filter propName within propRange
    F = false(size(obj.flatData.B,1),1);
    
    for k = 1:numel(obj.lineFiles)
        A = obj.loadLineK(k);
        cellIDs = A(:,obj.a('cellID'));
        P = obj.getPropAIndex((1:size(A,1))',A,propName,k);
        tf = (P>=propRange(1)) & (P<=propRange(2));
        
        z = accumarray(cellIDs,tf,[max(cellIDs) 1]);
        n = accumarray(cellIDs,1,[max(cellIDs) 1]);
        
        n(n==0) = eps;
        
        cellPassed = find(z==n);
        
        [fov,line] = obj.fovLineK(k);
        
        bIDs = obj.getBIndex(fov,line,cellPassed); %#ok<FNDSB>
        
        F(bIDs) = true;
    end
    
    end
    
    
    function F = filtersInUse(obj)
    %Retrieves a list of the filters and their ranges that are currently in
    %use
    F = struct;
    
    filterNames = fields(obj.cellFilterFlat);
    for k = 1:numel(filterNames)
        filterRange = obj.cellFilterFlat.(filterNames{k});
        if (numel(filterRange)>1 && any(isfinite(filterRange))) || (numel(filterRange)==1 && isfinite(filterRange) && filterRange~=0)
            F.(filterNames{k}) = filterRange;
        end
    end
    
    filterNames = fields(obj.cellFilterLine);
    for k = 1:numel(filterNames)
        filterRange = obj.cellFilterLine.(filterNames{k});
        if (numel(filterRange)>1 && any(isfinite(filterRange))) || (numel(filterRange)==1 && isfinite(filterRange) && filterRange~=0)
            F.(filterNames{k}) = filterRange;
        end
    end

    end

    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% FLAT: PARENT/CHIlD MAPPING
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    function obj = initParentList(obj)
    %OBJ = OBJ.GETPARENTLIST
    %Resets and runs obj.getChildren to retrieve a list of parents and
    %children.
    
    obj.parentList = {};
    obj.getChildren;
    
    end
    
    
    function [parentsID,childrenID] = getChildren(obj)
    %[parentsID,childrenID] = GETCHILDREN(obj).
    %Retrieves a list of parent indices (corresponds to obj.flatData.B
    %rows), and a cell with the indices of their respective children.
    
    if isempty(obj.parentList)
        [cID,pID] = obj.mkUniqueID;
        parents = unique(pID(pID~=0));
        children = cell(size(parents,1),1);
        
        for k = 1:size(parents,1);
            children{k} = cID(pID==parents(k))';
        end
        
        childrenID = children;
        parentsID = parents;
        
        obj.parentList = {};
        obj.parentList.parentsID = parentsID;
        obj.parentList.childrenID = childrenID;
    else
        parentsID = obj.parentList.parentsID;
        childrenID = obj.parentList.childrenID;
    end
    
    end
    
    
    
    function wholeCC = getWholeCC(obj)
    %wholeCC = GETWHOLECC(obj)    
    %Retrieves the indices in B of whole cell cycles.
        [parentsID,childrenID] = obj.getChildren;
        childrenID = cell2mat(reshape(childrenID,1,numel(childrenID)));
        childrenID = childrenID(:);
        wholeCC = childrenID(ismember(childrenID,parentsID));
    end
    
    function wholeCC = getWholeCCRobust(obj)
    %wholeCC = GETWHOLECCROBUST(obj)
    %Retrieves the indices of whole cell cycles after filtration has been
    %set such that only cells whose respective parents and children 
    %(in addition to themselves) passed the filters are kept. 
        [~,pID1] = obj.mkUniqueID;
        [pID2,cID2] = obj.getChildren;
        wholeCC0 = obj.getWholeCC;
        
        %Self check
        selfOK = obj.filterIDFlat(wholeCC0);
        wholeCC0 = wholeCC0(selfOK);
        
        %Parent check
        parentsID = pID1(wholeCC0);
        parentOK = obj.filterIDFlat(parentsID);
        
        wholeCC0 = wholeCC0(parentOK);
        
        %Child check
        childOK = false(numel(wholeCC0),1);
        for j = 1:numel(wholeCC0)
            childID = cID2{pID2==wholeCC0(j)};
            childOK(j) = all(obj.filterIDFlat(childID));
        end        
        
        wholeCC = wholeCC0(childOK);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% BUILDING UNIQUE CELL IDENTITIES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    function obj = initUniqueID(obj)
        cellID = (1:size(obj.flatData.B,1))';
        
        if min(obj.flatData.B(:,obj.b('fov')))==0
            addFov = 1;
        else
            addFov = 0;
        end
                
        matSize = max(obj.flatData.B(:,[obj.b('fov') obj.b('line') obj.b('cellID')]));
        matSize(1) = matSize(1) + addFov;
        C = zeros(matSize);
        cc = sub2ind(matSize,obj.flatData.B(:,obj.b('fov'))+addFov,obj.flatData.B(:,obj.b('line')),obj.flatData.B(:,obj.b('cellID')));
        C(cc) = cellID;
                
        obj.uniqueID = {};
        obj.uniqueID.ID = C;
        obj.uniqueID.fov0 = addFov;    
    end
    
    function [cellID,parentID] = mkUniqueID(obj)
        cellID = (1:size(obj.flatData.B,1))';
        
        isChild = obj.flatData.B(:,obj.b('parentID'));
        isChild = isChild~=0;
        
        pp = sub2ind(size(obj.uniqueID.ID),obj.flatData.B(isChild,obj.b('fov'))+obj.uniqueID.fov0,obj.flatData.B(isChild,obj.b('line')),...
            obj.flatData.B(isChild,obj.b('parentID')));
        
        parentID = 0*cellID;
        parentID(isChild) = obj.uniqueID.ID(pp);
    end
    
    
    
    function bID = getBIndex(obj,fov,line,lineID)
    %bID = obj.getBIndex(fov,line,lineID)
    %Get the index in B for a cell according to its field of view, line,
    %and cellID within the line. fov and line can either be constants, with
    %multiple lineIDs or 1-1 with lineID. Note bID is returned in the same
    %shape as lineID
    
    lineID0 = lineID;
    lineID = lineID(:);
    lineID = lineID(lineID~=0 & isfinite(lineID));
    
    if numel(fov)==1 && numel(line)==1 
        fov = fov*ones(numel(lineID),1);
        line = line*ones(numel(lineID),1);
    elseif (numel(fov)~=1 && numel(fov)~=numel(lineID)) || (numel(line)~=1 && numel(line)~=numel(lineID))
        warning('fov or line do not uniquely determine cell idenitities.')
        bID = 0*lineID(:);
        return
    end
        
    id = sub2ind(size(obj.uniqueID.ID),fov+obj.uniqueID.fov0,line,lineID);
    bID = obj.uniqueID.ID(id);
    
    bID0 = 0*lineID0;
    bID0(lineID0~=0 & isfinite(lineID0)) = bID;
    bID = bID0;
    
    end
    
    function bID = aID2bID(obj,K,aID,A)
    %bID = aID2bID(obj,K,aID,A)
    %Determines the index in B corresponding to the indices in A for line
    %K. K can be a vector but it must be the same size as aID.
    
    [fov,line] = obj.fovLineK(K);
    cellID = A(aID,obj.a('cellID'));
    
    bID = obj.getBIndex(fov,line,cellID);
    
    end
    
    function [K,aID] = bID2aID(obj,bID)
    %[k,aID] = bID2aID(obj,bID).
    %Determines the line and indices in matrix A corresponding to the
    %indices in bID
    fovLineCell = obj.flatData.B(bID,[obj.b('fov') obj.b('line') obj.b('cellID')]);
    matname = sprintf('%.2d_%.4d.mat',fovLineCell(:,[1 2])');
    matname = (reshape(matname,11,numel(matname)/11))';
    [~,K] = ismember(matname,obj.lineFiles);
        
    aID = cell(numel(bID),1);
    uniqueK = unique(K);
    for k = uniqueK'
        cellIDs = fovLineCell(K==k,3);
        A = obj.loadLineK(k);
        
        cell2aid = obj.cell2aid(A,k,false);
                
        aID(K==k) = cell2aid(cellIDs);
    end
    
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% MAP CONSECUTIVE A MATRIX INDICES 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function doaIDmap(obj)
    M = struct;
        for k = 1:numel(obj.lineFiles)
            A = obj.loadLineK(k);
            aID = (1:size(A,1))';
            
            [nxt,genUsed] = obj.getNextCyclePtaID(k,aID,A,'next');
            M(k).next = nxt;
            M(k).nextgenUsed = genUsed;
            
            [prev,genUsed] = obj.getNextCyclePtaID(k,aID,A,'prev');
            M(k).prev = prev;
            M(k).prevgenUsed = genUsed;
            
            M(k).last0 = obj.getEndCyclePtaID(0,aID,A,'last',0);
            M(k).last1 = obj.getEndCyclePtaID(0,aID,A,'last',1);
            
            M(k).first0 = obj.getEndCyclePtaID(0,aID,A,'first',0);
            M(k).first1 = obj.getEndCyclePtaID(0,aID,A,'first',1);
        end
        
        obj.aIDmap = M;
        obj.aIDmap(1).built = true;
        
    end
        
    %%% GET SPECIFIC CELL PROPERTIES
    function P = flat_cellProps(obj,prop)
        if ~iscell(prop), prop = {prop}; end
        
        P = {};
        
        for k = 1:numel(prop)
            P.(prop{k}) = obj.getPropBIndex(obj.filterIDFlat,prop{k});
        end        
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% CELL: RETRIEVE A SPECIFIC CELL
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    function cellData = getSingleCell(obj,bID)
        %Get the data from a single cell according to its index in B
        cellofInterest = obj.flatData.B(bID,:);
        cellData.B = cellofInterest;
        
        fov = cellofInterest(obj.b('fov'));
        line = cellofInterest(obj.b('line'));
        lineID = cellofInterest(obj.b('cellID'));
        
        lineInfo = load(fullfile(obj.dirList.lineDir,sprintf('fov%.2d',fov),sprintf('line%.4d',line),sprintf('%.2d_%.4d_INFO.mat',fov,line)));
        cellDataMAT = load(fullfile(obj.dirList.lineDir,sprintf('fov%.2d',fov),sprintf('line%.4d',line),[lineInfo.cellLine{lineID} '.mat']));
        
        cellData.cellData = cellDataMAT.cellData;
    end
    
    
    
    function bID = getBIndexCellData(obj,cellData)
        %Get the index in B from a cell according to its cellDataMP
        %structure
        fov = str2double(cellData.label(1:2));
        line = str2double(cellData.label(4:7));
        lineID = cellData.lineID;
        bID = obj.getBIndex(fov,line,lineID);
    end
        
    function cellIntensity = getCellIntensityVals(obj,cellData)
        %Get all intensity values corresponding to a given cell.
        %cellIntensity is a cell with elements corresponding to each frame
        %that the cell is present.
        fov = str2double(cellData.label(1:2));
        tc = 1;
        cellIntensity = {};
        
        for t = cellData.frames(1):cellData.frames(2)
            C = load(fullfile(obj.dirList.mat2Dir,sprintf('fov%.2d_t%.4d.mat',fov,t)));
            [im,bl] = obj.imList.getfilenamesFOVT(fov,t,1);
            I = imread(im);
            I = double(I);
            Iblank = imread(bl);
            Iblank = double(Iblank);
            I = I - Iblank;
            I(I<0) = 0;
            
            boundaryPts = C.CELLS.boundary{cellData.cellID(tc)};
            cellIntensity{tc} = getCellIntensityMP(boundaryPts,I); %#ok<AGROW>
            tc = tc + 1;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% GUIs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %CHECK UP ON SPECIFIC FRAMES
    function obj = doGUI(obj,fov,t)
    %obj = obj.doGUI(fov,t)
    %Bring up GUI to evaluate a particular field of view, starting at frame
    %t
        if isempty(obj.imList)
            obj.imList = imListMP(obj.dirList.dataDir);
        end
        evalGUIparent(obj.dirList,fov,t,obj.imList)
    end
    
    
       
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% LINE: RETREIVE PROPERTIES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%% CALCULATE INSTANTANEOUS GROWTH RATE
    function [alpha,aID,familyMap,newCellFlag] = alphaInst(obj,K,prop,nPts)
        %Calculates the instantaneous growth rate for the Kth line in
        %obj.lineFiles.
        [aID,familyMap,newCellFlag] = buildTreeFiltered(obj,K,1,1);
        A = load(fullfile(obj.dirList.flatDir,obj.lineFiles{K}));
        A = A.A;
                
        alpha = cell(1,numel(aID));
        for j = 1:numel(aID)
            if strcmp(prop,'V')
                PL = A(aID{j},obj.a('L'));
                Pw = A(aID{j},obj.a('w'));
                P = 4/3*pi*(Pw/2).^3 + pi*(Pw/2).^2.*(PL-Pw);
            else
                P = A(aID{j},obj.a(prop));
            end
            T = A(aID{j},obj.a('t'));
            
            newCell = logical(newCellFlag{j});
            dP = 0*P;
            dP(newCell) = 1;
            dP(1) = 0;
            dP = cumsum(dP);
            
            P = P.*2.^dP;
            P0 = P;
            
            %Use sliding window of nPts
            %Build sliding window index
            sID = repmat(1:numel(aID{j}),nPts,1);
            sID = sID + (-floor(nPts/2):nPts-floor(nPts/2)-1)'*ones(1,numel(aID{j}));
            sID(sID<1) = 1;
            sID(sID>numel(aID{j})) = numel(aID{j});
            
            T2 = T(sID);
            T2 = T2 - ones(nPts,1)*T2(1,:);
            T3 = mean(T2,2);
            
            
            %Estimate exponential growth by linearization
            P = P(sID);
            y = log(P);
            x = [T3(:) ones(nPts,1)];
            pp = x\y;
            
%             %Compare against exponential fit
%             pp2 = zeros(1,size(P,2));
%             for m = 1:size(P,2)
%                 f0 = fit(T2(:,m),P(:,m),'exp1');
%                 pp2(m) = f0.b;
%             end
            
            %Correct the calculation at the beginning and end of each
            %lineage (since fewer points should be used)
            for m = 1:floor(nPts/2)
                ym = P(nPts-m:end,m);
                xm = [T2(nPts-m:end,m) ones(size(ym,1),1)];
%                 f0 = fit(xm(:,1),ym,'exp1');
%                 pp2(m) = f0.b;
                ppm = xm\log(ym);
                pp(:,m) = ppm;
            end
            
            mc = nPts-1;
            for m = size(sID,2)-nPts+floor(nPts/2)+2:size(sID,2)
                ym = P(1:mc,m);
                xm = [T2(1:mc,m) ones(size(ym,1),1)];
%                 f0 = fit(xm(:,1),ym,'exp1');
%                 pp2(m) = f0.b;
                ppm = xm\log(ym);
                pp(:,m) = ppm;
                mc = mc - 1;
            end
            
            pp = pp(1,:);
            
            alpha{j} = [T(:) pp(:)];
        end        
    end
    
    function [P,aID,familyMap,newCellFlag] = line_cellProps(obj,K,prop)
    %[P,aID,familyMap,newCellFlag] = line_cellProps(obj,K,prop);
    %Retrieves a given property for the Kth lineage.    
        catEndPts = 0;
        groupMode = 'indep';
        useFilt = 1;
        
        [aID,familyMap,newCellFlag,A] = obj.buildTree2(K,groupMode,useFilt,catEndPts);
        
        if ~iscell(prop), prop = {prop}; end
        
        P = {};
        for p = 1:numel(prop)
            P.(prop{p}) = cell(1,numel(aID));
        end
        
        for j = 1:numel(aID)
            for p = 1:numel(prop)
                P.(prop{p}){j} = obj.getPropAIndex(aID{j},A,prop{p},K);
            end
        end
    end
    
    
    
    function P = line_cellPropsSingleCell(obj,K,prop,groupByCell,catEndPts,child2use)
    %P = OBJ.PROPLINESINGLECELL(k,prop,groupByCell,catEndPts,child2use)
    %Retrieves the properties specified by prop. If groupByCell is set to 1
    %(true), P is output as a cell array with each cell corresponding to a
    %whole cell cycle; otherwise it is a flattened vector of ordered cell
    %cycles. Set catEndPts to 1 to concatenate the last parent point and 
    %first daughter points to each single cell trajectory. If catEndPts and
    %groupByCell, child2use can be set to either 'mean' or 'rand' to either 
    %take the mean of two daughters or select a random daughter; if not
    %groupByCell, a random daughter is chosen.
    
    if ~iscell(prop), prop = {prop}; end

    useFilter = 1;
    [aID,A,genFactor] = obj.wholeCycleID(K,groupByCell,catEndPts,useFilter);
    
    P = {};
    for p = 1:numel(prop)
        if ~isempty(aID)
            P.(prop{p}) = cell(1,numel(aID));
        else
            P.(prop{p}) = [];
        end
    end
    
    if isempty(aID), return; end
    
    if groupByCell && catEndPts
        %Data is separated by cellIDs and end points have been
        %concatenated and a single cell may have two trajectories, one
        %corresponding to each daughter at the division point (last point)
        for j = 1:numel(aID)
            for p = 1:numel(prop)
                if ~obj.isExtensive(prop{p}) %~any(strcmp(prop{p},{'V','Vhemiellip','Vcyl','area','L','Itot'}))
                    %Extensive property: if parent or daughter values are used,
                    %need to normalize the property using genFactor. Otherwise
                    %just use factor of 1.
                    maxMin = [1 1];
                else
                    maxMin = [0.5 2];
                end
                                
                if catEndPts
                    %Parent and child data represented at the ends
                    switch child2use
                        case 'mean'
                            P1 = obj.getPropAIndex(aID{j}(:,1),A,prop{p},K).*max(maxMin(1),min(maxMin(2),genFactor{j}(:,1)));
                            P2 = obj.getPropAIndex(aID{j}(:,end),A,prop{p},K).*max(maxMin(1),min(maxMin(2),genFactor{j}(:,end)));
                            P0 = mean([P1 P2],2);
                        case 'rand'
                            r = randi(numel(size(aID{j},2)),1,1);
                            P0 = obj.getPropAIndex(aID{j}(:,r),A,prop{p},K).*max(maxMin(1),min(maxMin(2),genFactor{j}(:,r)));
                    end
                else
                    P0 = obj.getPropAIndex(aID{j},A,prop{p},K);
                end
                P.(prop{p}){j} = P0;
            end
        end
        
    elseif groupByCell && ~catEndPts
        %No daughter averages to calculate, can treat as a flat vector
        nC = cellfun(@numel,aID);
        aID0 = cell2mat(aID(:));
        
        for p = 1:numel(prop)
            %N.B. no extensivity check is required since no indices cross generations 
            P0 = obj.getPropAIndex(aID0,A,prop{p},K);
            P.(prop{p}) = mat2cell(P0,nC,1);
        end
    else
        %Data is contained in a flattened vector
        for p = 1:numel(prop)
            if ~obj.isExtensive(prop{p}) %~any(strcmp(prop{p},{'V','Vcyl','area','L','Itot'}))
                %Extensive property: if parent or daughter values are used,
                %need to normalize the property using genFactor. Otherwise
                %just use factor of 1.
                genFactorP = 1;%ones(numel(genFactor),1);
            else
                genFactorP = genFactor;
            end
            
            P.(prop{p}) = obj.getPropAIndex(aID,A,prop{p},K).*genFactorP;
        end
    end
    
    end
    
    
    
    
    function [P,familyData,binValue] = line_cellPropsGenerationT(obj,K,T,firstGen,lastGen,prop,binProp,leftRight)
    %[P,familyData,binID] = obj.line_cellPropsGenerationT(K,T,firstGen,lastGen,prop,binProp,leftRight)
    %Get properties specified by prop for the Kth line, for the firstGen
    %through lastGen from time point T. P is a structure with fields 
    %corresponding to prop. If a field refers to an instantaneous property 
    %(from A matrices), it retains the same size and shape as the list of 
    %indices (aID) output by getGenerationsT. If the field refers to a flat
    %property (from B matrix), it has the same size and shape as familyMap,
    %with NaNs where no cell in the lineage was tracked (or if it did not 
    %pass the filter). binValue contains the value of the
    %properties in binProp at the time nearest to T for each lineage.
    %leftRight determines whether to take the value to the 'left', 'right',
    %or 'any' side of time T.
    
    %Get a list of the relevant generations from time T
    useFilt = 1;
    catEndPts = 1;
    [aID,familyMap,newCellFlag,gen0IDs,genMap,A] = obj.getGenerationsT(K,T,firstGen,lastGen,useFilt,catEndPts);
    
    familyData = {};
    familyData.aID = aID;
    familyData.familyMap = familyMap;
    familyData.newCellFlag = newCellFlag;
    familyData.genMap = genMap;
    familyData.gen0IDs = gen0IDs;
        
    if ~iscell(prop); prop = {prop}; end    
    if ~iscell(binProp); binProp = {binProp}; end
    
    %Initialize property structure
    P = {};
    for p = 1:numel(prop)
        P.(prop{p}) = cell(size(aID));
    end    
    
    %Split properties depending on matrix A or B
    [propA,propB] = obj.propParse(prop);
    [binPropA,binPropB] = obj.propParse(binProp);
    
    if ~isempty(propB) || ~isempty(binPropB)
        [fov,line] = obj.fovLineK(K);
    end
    
    %Save properties
    for j = 1:numel(aID)
        %Properties from A
        for p = 1:numel(propA)
            P.(propA{p}){j} = cell(size(aID{j}));
            for m = 1:numel(aID{j})
                P.(propA{p}){j}{m} = obj.getPropAIndex(aID{j}{m},A,propA{p});
            end
        end
        
        %Properties from B
        cellIDs = unique(familyMap{j}(:));
        cellIDs = cellIDs(cellIDs~=0);    
        for p = 1:numel(propB)
            P.(propB{p}){j} = nan(size(familyMap{j}));
            bID = obj.getBIndex(fov,line,cellIDs);
            for c = 1:numel(cellIDs)
                if ~bID(c)
                    %Index not found in B
                    continue
                else
                    P.(propB{p}){j}(familyMap{j}==cellIDs(c)) = obj.getPropBIndex(bID(c),propB{p});
                end
            end
        end
    end
    
    %Values for binning by the zeroth generation
    [aID0,cc0] = obj.getAIndexT(T,leftRight,A,gen0IDs);
        
    binValue = {};
    binValue.cellCycleT = cc0;
    
    %Find values in A and B matrices
    for p = 1:numel(binPropA)
        binValue.(binPropA{p}) = obj.getPropAIndex(aID0,A,binPropA{p});
    end
    for p = 1:numel(binPropB)
        binValue.(binPropB{p}) = zeros(numel(gen0IDs),1);
        for c = 1:numel(gen0IDs)
            binValue.(binPropB{p})(c) = obj.getPropBIndex(obj.getBIndex(fov,line,gen0IDs(c)),binPropB{p});
        end
    end
    
    end
    
    
    
    
    function [P,familyData,binValue] = line_cellPropsGenerationT2(obj,K,T,firstGen,lastGen,prop,binProp,leftRight)
    %[P,familyData,binID] = obj.line_cellPropsGenerationT2(K,T,firstGen,lastGen,prop,binProp,leftRight)
    %Get properties specified by prop for the Kth line, for the firstGen
    %through lastGen from time point T. P is a structure with fields 
    %corresponding to prop. If a field refers to an instantaneous property 
    %(from A matrices), it retains the same size and shape as the list of 
    %indices (aID) output by getGenerationsT2. If the field refers to a flat
    %property (from B matrix), it has the same size and shape as familyMap,
    %with NaNs where no cell in the lineage was tracked (or if it did not 
    %pass the filter). binValue contains the value of the
    %properties in binProp at the time nearest to T for each lineage.
    %leftRight determines whether to take the value to the 'left', 'right',
    %or 'any' side of time T.
    
    %Get a list of the relevant generations from time T
    useFilt = 1;
    catEndPts = 0;
    [aID,familyMap,newCellFlag,gen0IDs,genMap,genaID,A] = obj.getGenerationsT2(K,T,firstGen,lastGen,useFilt,catEndPts);
    
    familyData = {};
    familyData.aID = aID;
    familyData.familyMap = familyMap;
    familyData.newCellFlag = newCellFlag;
    familyData.genMap = genMap;
    familyData.gen0IDs = gen0IDs;
    familyData.genaID = genaID;
        
    if ~iscell(prop); prop = {prop}; end    
    if ~iscell(binProp); binProp = {binProp}; end
    
    %Initialize property structure
    P = {};
    for p = 1:numel(prop)
        P.(prop{p}) = cell(size(aID));
    end    
    
    %Split properties depending on matrix A or B
    [propA,propB] = obj.propParse(prop);
    [binPropA,binPropB] = obj.propParse(binProp);
    
    if ~isempty(propB) || ~isempty(binPropB)
        [fov,line] = obj.fovLineK(K);
    end
    
    %Save properties
    for j = 1:numel(aID)
        %Properties from A
        for p = 1:numel(propA)
            P.(propA{p}){j} = obj.getPropAIndex(aID{j},A,propA{p});
        end
        
        %Properties from B
        for p = 1:numel(propB)
            bID = obj.getBIndex(fov,line,familyMap{j});
            P.(propB{p}){j} = obj.getPropBIndex(bID,propB{p});
        end
    end
    
    %Values for binning by the zeroth generation
    [aID0,cc0] = obj.getAIndexT(T,leftRight,A,gen0IDs);
        
    binValue = {};
    binValue.cellCycleT = cc0;
    
    %Find values in A and B matrices
    for p = 1:numel(binPropA)
        binValue.(binPropA{p}) = obj.getPropAIndex(aID0,A,binPropA{p});
    end
    
    for p = 1:numel(binPropB)
        binValue.(binPropB{p}) = obj.getPropBIndex(obj.getBIndex(fov,line,gen0IDs),binPropB{p});
    end
    
    end
    
    
    
    
      
    function [aID,cellcycleT] = getAIndexT(obj,T,leftRight,A,cellIDs)
    %[aID,cellcycleT] = obj.getAIndexT(T,leftRight,A,cellIDs)
    %Get the index in A for each cell in cellIDs that is closest to the
    %timepointT, approaching either from the 'left', 'right', or 'any'
    %direction. cellcycleT contains the time since the cell's birth (first
    %time point observed) in the first column and the fraction of its
    %(observed) cell cycle that has passed by that point in the second 
    %column.
    
    aID = zeros(numel(cellIDs),1);
    cellcycleT = zeros(numel(cellIDs),2);
    for j = 1:numel(cellIDs)
        aIDj = find(A(:,obj.a('cellID'))==cellIDs(j));
        tj = A(aIDj,obj.a('t'));
        tj = tj - T;
        switch leftRight
            case 'left'
                tj(tj>0) = Inf;
            case 'right'
                tj(tj<0) = Inf;                
        end
        tj = abs(tj);
        aID(j) = aIDj(find(tj==min(tj),1,'first'));
        cellcycleT(j,:) = [A(aID(j),obj.a('t'))-A(aIDj(1),obj.a('t')) find(tj==min(tj),1,'first')./length(tj)];
    end
    
    end        
       
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% LINE: LINEAGE RECONSTRUCTIONS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    function [aID,familyMap,newCellFlag] = buildTreeFiltered(obj,K,indep,catEndPts)
        %Filter using the same filter that is currently set for flat
        %analysis. Cell lines are broken if a cell should be filtered out,
        %but the first and last points of its cell cycle are kept in each
        %Build unfiltered tree
        [aID,familyMap,newCellFlag] = obj.buildTree(K,indep);
        
        fov = str2double(obj.lineFiles{K}(1:2));
        line = str2double(obj.lineFiles{K}(4:7));
        
        %Filter using the cells status in the flat filter index
        k = 1;
        while k <= numel(aID)
            cellIDs = familyMap{k};
            rmCell = zeros(numel(cellIDs),1);
            for j = 1:numel(cellIDs)
                bID = obj.getBIndex(fov,line,cellIDs(j));
                rmCell(j) = not(obj.filterIDFlat(bID));
            end
            rmCell = find(rmCell);
            for j = 1:numel(rmCell)
                if rmCell(j)==1 && numel(familyMap{k})>1 
                    %Remove the first cell in the series
                    f1 = find(newCellFlag{k}==cellIDs(rmCell(j)+1));
                    aID{k} = aID{k}(f1:end);
                    newCellFlag{k} = newCellFlag{k}(f1:end);
                    familyMap{k} = newCellFlag{k}(newCellFlag{k}~=0);
                
                elseif cellIDs(rmCell(j))==familyMap{k}(end);
                    %Remove the last cell in the series
                    f2 = find(newCellFlag{k}==cellIDs(rmCell(j)))-1;
                    aID{k} = aID{k}(1:f2);
                    newCellFlag{k} = newCellFlag{k}(1:f2);
                    familyMap{k} = newCellFlag{k}(newCellFlag{k}~=0);
                    
                elseif rmCell(j)==1 && numel(familyMap{k})==1
                    %Remove the entire lineage
                    aID(k) = [];
                    newCellFlag(k) = [];
                    familyMap(k) = [];
                    break
                else
                    %Break up the lineage: remove the relevant cell and add
                    %the remainder of the line to a new cell
                    f3 = find(newCellFlag{k}==cellIDs(rmCell(j)))-1;
                    f4 = find(newCellFlag{k}==cellIDs(rmCell(j)+1));
                    
                    %Store the rest of the line at the end of the cell
                    naID = numel(aID);
                    aID{naID+1} = aID{k}(f4:end);
                    newCellFlag{naID+1} = newCellFlag{k}(f4:end);
                    familyMap{naID+1} = newCellFlag{naID+1}(newCellFlag{naID+1}~=0);
                    
                    %Truncate the current index
                    aID{k} = aID{k}(1:f3);
                    newCellFlag{k} = newCellFlag{k}(1:f3);
                    familyMap{k} = newCellFlag{k}(newCellFlag{k}~=0);
                    
                    %Break from current loop: algorithm will deal with the
                    %rest of lineage through the while loop (k)
                    break
                end
            end
            k = k + 1;
        end
        
        anyLeft = logical(cellfun(@numel,aID));
        aID = aID(anyLeft);
        newCellFlag = newCellFlag(anyLeft);
        familyMap = familyMap(anyLeft);
        
        if catEndPts
            %Concatenate the endpoints 
            A = load(fullfile(obj.dirList.flatDir,obj.lineFiles{K}));
            A = A.A;
            
            [aID,familyMap,newCellFlag] = buildTreeEndPtCat(obj,A,aID,familyMap,newCellFlag);            
        end
        
    end
    
    function [aID,familyMap,newCellFlag] = buildTreeEndPtCat(obj,A,aID,familyMap,newCellFlag)
        %Concatenates the first and last points of each family line with the last and first points of the
        %corresponding parent and daughter(s). If there are two daughters,
        %one is chosen at random. Concatenation is performed on all three
        %aID, familyMap, and newCellFlag.
        
        for k = 1:numel(aID)                
            hasParent = A(aID{k}(1),obj.a('parentID'));
            if hasParent~=0
                aIDparent = find(A(:,obj.a('cellID'))==hasParent,1,'last');
                aID{k} = [aIDparent; aID{k}];
                newCellFlag{k} = [hasParent; newCellFlag{k}];
                familyMap{k} = [hasParent; familyMap{k}];
            end
            hasChild = A(:,obj.a('parentID'))==A(aID{k}(end),obj.a('cellID'));
            hasChild = A(hasChild,obj.a('cellID'));
            if numel(hasChild)>0
                rC = randi(numel(hasChild),1,1);
                hasChild = hasChild(rC);
                aIDchild = find(A(:,obj.a('cellID'))==hasChild,1,'first');
                aID{k} = [aID{k}; aIDchild];
                newCellFlag{k} = [newCellFlag{k}; hasChild];
                familyMap{k} = [familyMap{k}; hasChild];
            end
        end
    end
    
    
    
    function [aID,familyMap,newCellFlag] = buildTree(obj,k,indep)
        %[aID,familyMap,newCellFlag] = buildTree(obj,k,indep)
        %Reconstructs cell lines as vectors containing the (chronological) 
        %indices in the aID-th line file listed in obj.lineFiles. Set
        %indep to 0 or 1 to list lines in their entirety, including common
        %parents (0) or independent lines (1).
        
        A = load(fullfile(obj.dirList.flatDir,obj.lineFiles{k}));
        A = A.A;
        if indep
            %Save only the independent lines, starting with the longest lines 
            familyMap = {};
            aID = {};
            newCellFlag = {};
            c = 1;
            cellIDs = unique(A(:,obj.a('cellID')));
            cell2written = zeros(max(cellIDs),1);
            cell2written(cellIDs) = 1:numel(cellIDs);
            cellWritten = false(numel(cellIDs),1);
            while any(~cellWritten)
                %Look for the line of maximum length
                cellsLeft = cellIDs(~cellWritten)';
                line = cell(numel(cellsLeft),1);
                for j = 1:numel(cellsLeft)
                    linej = obj.getLineage(A(:,[obj.a('cellID') obj.a('parentID')]),cellsLeft(j),cellIDs(cellWritten));
                    alreadyWritten = cellWritten(cell2written(linej));
                    if sum(alreadyWritten)==0
                        line{j} = linej;
                    else
                        line{j} = linej(1:find(alreadyWritten,1)-1);
                    end
                end
                
                maxGen = cellfun(@numel,line);
                [~,maxID] = max(maxGen);
                
                familyMap{c} = line{maxID}(end:-1:1);
                aIDc = [];
                for j = 1:numel(familyMap{c})
                    aIDc = [aIDc; find(A(:,obj.a('cellID'))==familyMap{c}(j))];
                end
                aID{c} = aIDc(:);
                newCellFlag{c} = A(aIDc,obj.a('cellID')).*[true; diff(A(aIDc,obj.a('cellID')))~=0];
                
                cellWritten(cell2written(line{maxID})) = true;
                c = c + 1;
            end            
        else
            %Get all lines, even redundant (common) cell indices
            cellIDs = unique(A(:,obj.a('cellID')));
            familyMap = cell(1,numel(cellIDs));
            aID = cell(1,numel(cellIDs));
            newCellFlag = cell(1,numel(cellIDs));
            for j = 1:numel(cellIDs)
                familyMap{j} = obj.getLineage(A(:,[find(obj.a('cellID')) find(obj.a('parentID'))]),cellIDs(j),Inf);
                aIDc = [];
                for c = 1:numel(familyMap{j})
                    aIDc = [aIDc; find(A(:,obj.a('cellID'))==familyMap{j}(c))];
                end
                aID{j} = aIDc(:);
                newCellFlag{j} = A(aIDc(:),obj.a('cellID')).*[true; diff(A(aIDc(:),obj.a('cellID')))~=0];
            end
        end        
    end
    
    function [aID,familyMap,newCellFlag,A] = buildTree2(obj,K,groupMode,useFilt,catEndPts)
    %[aID,familyMap,newCellFlag,A] = buildTree2(obj,K,groupMode,useFilt,catEndPts)
    %Builds the family tree for the K-th line file. Set groupMode to either
    %'indep' or 'family' to generate either independent or lineage-grouped
    %lines. If lineage-grouping is used, note that the outputs will be
    %matrices with one column per unique daughter. Set useFilt to ignore
    %any cells that would be removed by the current filter. Set catEndPts
    %to concatenate the first parent and last daughter cell cycle points. 
    %Note that concatenated parents and daughters are added to the new cell
    %flag but not the family map.
    
    A = obj.loadLineK(K);
    
    [aIDk,AcellID] = obj.cell2aid(A,K,useFilt);
    
    %Find the first cell of all possible lines
    firstGen = AcellID(:,2)==0 & AcellID(:,1)~=0;
    firstGen = unique(AcellID(firstGen,1));

    nMax = numel(firstGen);
    familyMap = cell(1,nMax);

    nLines = zeros(1,nMax);
    for k = 1:nMax
        lineK = obj.getLineageForward(AcellID,firstGen(k),nMax);
        lineEnd = find(sum(lineK,2)==0,1,'first');
        if ~isempty(lineEnd)
            lineK = lineK(1:lineEnd-1,:);
        end
        lineK = sortrows(lineK');
        familyMap{k} = lineK';

        nLines(k) = size(familyMap{k},2);
    end
    
    switch groupMode
        case 'indep'
            %Separates into independent lines
            familyMapI = cell(1,sum(nLines));
            c = 1;
            for k = 1:nMax
                fmk = familyMap{k};
                dfmk = [ones(size(fmk,1),1) diff(fmk,1,2)];
                dfmk = logical(dfmk) & fmk~=0;

                fmk = num2cell(fmk);
                fmk(~dfmk) = {[]};
                fmk = mat2cell(fmk,size(fmk,1),ones(1,nLines(k)));
                fmk = num2cell(fmk);
                fmk = cellfun(@(x) cell2mat(x{1}),fmk,'uniformoutput',0);

                familyMapI(c:c+nLines(k)-1) = fmk;
                c = c + nLines(k);
            end
            familyMap = familyMapI;

            %Build the corresponding indices
            aID = cell(size(familyMap));
            newCellFlag = cell(size(familyMap));
            for k = 1:numel(familyMap)
                aid = obj.familyMap2aID(familyMap{k},aIDk);
                aID{k} = aid;
                newCellFlag{k} = logical([1; diff(AcellID(aid,1))]).*AcellID(aid,1);
            end

        case 'family'
            %Keeps family mapping and indices grouped according to family
            aID = cell(size(familyMap));
            newCellFlag = cell(size(familyMap));
            for k = 1:numel(familyMap)
                aID{k} = obj.familyMap2aID(familyMap{k},aIDk);

                ncf = aID{k};
                ncf(ncf~=0) = AcellID(ncf(ncf~=0),1);
                newCellFlag{k} = ncf.*[ones(1,size(ncf,2)); diff(ncf)];  
            end
    end
    
    %Concatenate end points if necessary
    if catEndPts
        [aID,newCellFlag] = obj.catEndPts(K,A,aID,newCellFlag);
    end

    end
   
    
    
    
    
    
    
    function [aID,newCellFlag,genaID] = catEndPts(obj,K,A,aID,newCellFlag,genaID)
    %[aID,newCellFlag,genaID] = catEndPts(obj,K,A,aID,newCellFlag,genaID)
    %Concatenates end points onto aID, newCellFlag, and genaID (if supplied).
    
    flag = false;
    if nargin<6
        %If no generational aID is supplied, just create a dummy of the
        %correct size for manipulation.
        warning('No genaID was provided.')
        flag = true;
        genaID = aID;
    end
    
    for k = 1:numel(aID)
        %Concatenate indices
        aidk = aID{k};
        gk = genaID{k};
        
        getMyParent = aidk(1,:);
        parentaID = obj.getEndCyclePtaID(K,getMyParent,A,'first',1);
        parentaID(isnan(parentaID)) = [];

        getMyDaughter = [logical(aidk); zeros(1,size(aidk,2))];
        getMyDaughter = diff(getMyDaughter,1,1);
        getMyDaughter0 = find(getMyDaughter==-1);
        getMyDaughter = aidk(getMyDaughter0);

        daughteraID = obj.getEndCyclePtaID(K,getMyDaughter,A,'last',1);
        daughteraID(isnan(daughteraID)) = 0;
        
        getMyDaughter0 = getMyDaughter0(:);
        if ~all(daughteraID==0)
            aidk = [aidk; zeros(1,size(aidk,2))]; %#ok<AGROW>
            aidk(getMyDaughter0+(1:numel(getMyDaughter0))') = daughteraID;
            gk = [gk; zeros(1,size(gk,2))]; %#ok<AGROW>
            gk(getMyDaughter0+(1:numel(getMyDaughter0))') = genaID{k}(getMyDaughter0)+1;
        end
        
        parentaID = reshape(parentaID,numel(parentaID),1);
        
        aID{k} = [parentaID'; aidk];
        if ~flag
            genaID{k} = [gk(1,:)-1; gk];
        end
        
        %Write parent and daughters to newCellFlag
        ncfk = [newCellFlag{k}; zeros(1,size(aidk,2))];
        parentncf = A(parentaID,obj.a('cellID'));
        daughterncf = zeros(1,numel(daughteraID));
        daughterncf(daughteraID~=0) = A(daughteraID(daughteraID~=0),obj.a('cellID'));

        ncfk(getMyDaughter0+(1:numel(getMyDaughter0))') = daughterncf;

        newCellFlag{k} = [parentncf'; ncfk];            
    end
    
    end
    
    
    
    
    function P = singleCC_cellProps(obj,props,groupByCell,catEndPts)
    %P = singleCC_cellProps(obj,prop,groupByCell,catEndPts) 
    %Using dataMP2.line_cellPropsSingleCell, checks all cell lines to build 
    %the structure P with fields props, each of which is a cell array 
    %containing cell cycle values, where each cell in the array corresponds
    %to a single cell cycle.
    
    if ~iscell(props); props = {props}; end
    
    if nargin < 3
        groupByCell = true; 
        catEndPts = false;
    elseif nargin < 4
        catEndPts = false;
    end
    
    %Initialize cell array
    P = {};
    for p = 1:numel(props)
        if groupByCell
            P.(props{p}) = cell(1e4,1);
        else
            %Check the output size
            [~,SZ] = obj.shapeConstant(props{p},'A',[3 1]);
            P.(props{p}) = zeros(5e4,SZ(2));
        end
    end
       
    %Retrieve properties for all cell lines
    c = 1;
    for k = 1:numel(obj.lineFiles)
        Pk = obj.line_cellPropsSingleCell(k,props,groupByCell,catEndPts,'mean');
        nk = numel(Pk.(props{1}));
        
        %Store values
        for p = 1:numel(props)
            P.(props{p})(c:c+nk-1,:) = Pk.(props{p});
        end
        c = c + nk;
    end
    
    %Truncate where necessary
    for p = 1:numel(props)
        P.(props{p}) = P.(props{p})(1:c-1,:);
    end
    
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    function [aID,A,genFactor] = wholeCycleID(obj,K,groupByCell,catEndPts,useFilter)
    %[aID,A,genFactor] = wholeCycleID(obj,K,groupByCell,catEndPts,useFilter)
    %Retrieves the indices in the Kth lineFiles A matrix corresponding 
    %to whole cell cyles. If groupByCell is set to 1 (true), aID is an
    %array with cells corresponding to single cell cycles; otherwise aID is
    %a vector. Set catEndPts to 1 (true) to concatenate the parent and
    %daughter indices, note that if groupByCell is false only one daughter
    %is concatenated (chosen at random); otherwise, aID contains two 
    %columns. Set useFilter to true to filter cells based on the current 
    %flat filter in place. genFactor is the same shape as aID, and contains
    %the factor to use to halve or double extensive property values due to
    %the concatenation of daughter and parent indices.
    
    %Load data file
    A = load(fullfile(obj.dirList.flatDir,obj.lineFiles{K}));
    A = A.A;
    
    %Find the cells belonging to whole cell cycles, i.e. they have both a
    %parent and a daughter
    hasParent0 = A(:,obj.a('parentID'));
    hasParent = A(hasParent0~=0,obj.a('cellID'));
    hasParent = unique(hasParent);
    isParent = A(:,obj.a('cellID'));
    isParent = unique(isParent);
    isParent = isParent(ismember(isParent,hasParent0));
    
    wholeCC = hasParent(ismember(hasParent,isParent));
    
    fov = str2double(obj.lineFiles{K}(1:2));
    line = str2double(obj.lineFiles{K}(4:7));
    
    bID = obj.getBIndex(fov,line,wholeCC);
    
    %Filter cells if necessary
    if useFilter
        wholeCC = wholeCC(obj.filterIDFlat(bID));
    end
    
    if isempty(wholeCC)
        aID = [];
        genFactor = [];
        return
    end
    
    if groupByCell && catEndPts
        %Keep cell cycles separated, concatenate end points
        aID = cell(1,numel(wholeCC));
        genFactor = cell(1,numel(wholeCC));
        for j = 1:numel(wholeCC)
            f1 = find(A(:,obj.a('cellID'))==wholeCC(j));
            g1 = ones(numel(f1),1);
            
            if catEndPts
                parentID = A(A(:,obj.a('cellID'))==wholeCC(j),obj.a('parentID'));
                parentID = parentID(1);
                f0 = find(A(:,obj.a('cellID'))==parentID,1,'last');
                daughterID = A(A(:,obj.a('parentID'))==wholeCC(j),obj.a('cellID'));
                daughterID = [daughterID(1) daughterID(end)];
                if diff(daughterID)==0
                    %Just one daughter
                    f2 = find(A(:,obj.a('cellID'))==daughterID(1),1,'first');
                    aID{j} = [f0; f1(:); f2];
                    genFactor{j} = [0.5; g1(:); 2];
                else
                    %Two daughters
                    f21 = find(A(:,obj.a('cellID'))==daughterID(1),1,'first');
                    f22 = find(A(:,obj.a('cellID'))==daughterID(2),1,'first');
                    aID{j} = [[f0; f1(:); f21] [f0; f1(:); f22]];
                    genFactor{j} = [[0.5; g1(:); 2] [0.5; g1(:); 2]];
                end
            else
                aID{j} = f1(:);
                genFactor{j} = g1(:);
            end
        end
        
    elseif groupByCell && ~catEndPts
        %If parent and daughter values are ignored, it is faster to treat
        %all indices as a vector and separate out into cells at the end
        aID0 = find(ismember(A(:,obj.a('cellID')),wholeCC));
        
        newCell = A(aID0,obj.a('cellID'));
        newCell = [diff(newCell)~=0; true];
        newCell = find(newCell);
        newCell = [newCell(1); diff(newCell)];
        
        aID = mat2cell(aID0,newCell,1);
        genFactor = mat2cell(ones(numel(aID0),1),newCell,1);
        
    else
        %Flatten all indices to a single vector
        aID0 = find(ismember(A(:,obj.a('cellID')),wholeCC));
                
        if catEndPts
            newCell = logical([1; diff(A(aID0,obj.a('cellID')))]);
            parentaID = obj.getEndCyclePtaID(K,aID0(newCell),A,'first',1);
            daughteraID = obj.getEndCyclePtaID(K,aID0(newCell),A,'last',1);
            
            aID = zeros(numel(aID0)+numel(parentaID)+numel(daughteraID),1);
            genFactor = ones(numel(aID),1);
            
            aIDIndex = 1:numel(aID0);
            addMe1 = cumsum(newCell);
            addMe2 = [0; cumsum(newCell(2:end))];
            aIDIndex = aIDIndex(:) + addMe1 + addMe2;
            
            %Store whole cell cycles
            aID(aIDIndex) = aID0;
            
            %Store parents
            aIDIndex = logical([diff(aID)>0; 0]) & aID==0;
            aID(aIDIndex) = parentaID;
            genFactor(aIDIndex) = 1/2;
            
            %Store daughters
            aIDIndex = aID==0;
            aID(aIDIndex) = daughteraID;
            genFactor(aIDIndex) = 2;
        else
            %No need to concatenate points: output a single vector
            aID = aID0;
            genFactor = ones(numel(aID),1); %Since no points concatenated, only a single geenration present            
        end
    end
    
    end    
    
    
    function P = getPropAIndex(obj,aID,A,prop,K)
    %P = obj.getPropAIndex(aID,A,prop,K)
    %Retrieves the property prop from matrix A according to indices aID.
    %Note that P will have the same shape as aID.
    
    if isempty(aID)
        P = [];
        return
    end
    
    %Initialize matrices to store values
    aID0 = aID;
    PaID = nan(size(aID0));
    
    aID = aID(aID~=0 & isfinite(aID));
    
    if regexp(prop,'^raw\w+')
        %Raw values obtained from segmentation, i.e. no smoothing applied,
        %or a raw calculation without smoothing, i.e. V and SA
        baseProp = regexp(prop,'raw(\w+)','tokens');
        baseProp = baseProp{1}{1};
        
        if any(strcmp(baseProp,{'V','Vhemisph','Vhemiellip','Vcyl','SA','SAhemisph','SAcyl','SAhemiellip'}))
            %Raw volume or surface area calculation
            [baseProp,mthd] = obj.parseMethod(prop);
            
            %Only need area for ellipsoid geometry
            if ~strcmp(mthd,'hemiellip')
                area = [];
            else
                area = obj.getPropAIndex(aID(:),A,'area',K);
            end
            
            L = obj.getPropAIndex(aID(:),A,'L',K);
            w = obj.getPropAIndex(aID(:),A,'w',K);
                        
            if strcmp(baseProp(1),'V')
                P = obj.volume(L,w,area,mthd);
            else
                P = obj.surfacearea(L,w,area,mthd);
            end
            
        elseif strcmp(baseProp,'w')
            %Raw width calculation from area
            area = obj.getPropAIndex(aID(:),A,'area',K);
            L = obj.getPropAIndex(aID(:),A,'L',K);
            
            P = obj.width(area,L);
            
        elseif strcmp(baseProp,'L');
            %Raw length from segmentation
            P = A(aID,obj.a('L')) + 0.1067;
            
        elseif isKey(obj.a,baseProp)
            %Retrieve segmented value
            P = A(aID,obj.a(baseProp));
            
        end
        
    elseif any(strcmp(prop,{'area','L','w','V','Vhemisph','Vhemiellip','Vcyl',...
            'SA','SAhemisph','SAhemiellip','SAcyl'})) && obj.smoothingMethod.dosmooth
        %Measurements that might have a smoothing filter applied that
        %underlies other calculations
        %Parse the prop string for potential methods
        [baseProp,mthd] = obj.parseMethod(prop);
        
        if isfield(obj.smoothingMethod.A,baseProp)
            if ~isempty(obj.smoothingMethod.A.(baseProp))
                %Retrieve the smoothed property
                %Format required string to retrieve property from
                %getPropAIndex
                str = sprintf(obj.smoothingMethod.A.(baseProp),['raw' baseProp mthd]);
                P = obj.getPropAIndex(aID(:),A,str,K);
            else
                %No smoothing required
                P = obj.getPropAIndex(aID(:),A,['raw' baseProp mthd],K);
            end
        else
            %No smoothing required
            P = obj.getPropAIndex(aID(:),A,['raw' baseProp mthd],K);
        end
                
    elseif strcmp(prop(1),'V') 
        L = obj.getPropAIndex(aID(:),A,'L',K);
        wArea = obj.getPropAIndex(aID(:),A,'w',K);
        w = wArea;
        
        area = A(aID,obj.a('area'));
        P = obj.volume(L,w,area,prop(2:end));
    
    elseif strcmp(prop(1),'C')
        %Use method specified for volume calculation, with the default
        %(i.e. if method is left empty) being cyl + hemisph caps
        if numel(prop)==1
            V = obj.getPropAIndex(aID(:),A,'V',K);
        else
            V = obj.getPropAIndex(aID(:),A,['V' prop(2:end)],K);
        end
        
        P = A(aID,obj.a('Itot'))./V;
        
    elseif strcmp(prop,'wMean');
        wRaw = obj.getPropAIndex(aID(:),A,'w',K);
        wCC = obj.getPropAIndex(aID(:),A,'cellCycw',K);
        P = (wRaw + wCC)./2;
        
    elseif strcmp(prop,'wareadbexpsm')
        %Obtain width from the double exponential smoothed area and length
        %(smoothed from the beginning of the cell cycle)
        area = obj.getPropAIndex(aID,A,'dbexpsmooth2area',K);
        L = obj.getPropAIndex(aID,A,'dbexpsmooth2L',K);
        
        P = obj.width(area,L);
        
    elseif strcmp(prop,'wRaw')
        P = A(aID,obj.a('w')) + 0.1067;
    
    elseif strcmp(prop(1),'w')
        %Obtain width from the segmented area under spherocylindrical
        %assumption
        area = A(aID,obj.a('area'));
        L = A(aID,obj.a('L'));
        L = L + 0.1067;
        
        P = 2*(L-sqrt(L.^2-(4-pi)*area))./(4-pi);
        
    elseif strcmp(prop(1),'L')
        %Add pixel offset to length
        P = A(aID,obj.a('L')) + 0.1067;
    
    elseif any(regexp(prop,'^cellCyc\w+$'))
        %The cell cycle mean
        baseProp = regexp(prop,'^cellCyc(\w+)$','tokens');
        baseProp = baseProp{1}{1};
                
        %Get the full cell cycle indices for each aID
        cell2aid = obj.cell2aid(A,K,false);
        
        %Find unique cellIDs corresponding to aID
        cellIDs = unique(A(aID(:),obj.a('cellID')));
        
        %Bin each aID to its cell
        cellBins = [cellIDs; cellIDs(end)+1];
        [~,~,binID] = histcounts(A(aID(:),obj.a('cellID')),cellBins);
        
        %Get all relevant aID in ascending order of cellIDs
        fullaID = cell2aid(cellIDs);
        nID = cellfun(@numel,fullaID(:)); %Alternative way to obtain nID for indexing, slower but difference is negligible (3%)        
        fullaID = cell2mat(fullaID(:));
        
        %Build the index for accumarray to act on
        z = 0*fullaID;
        z(1) = 1;
        z(cumsum(nID(1:end-1))+1) = 1;
        z = cumsum(z);
                
        Pfull = obj.getPropAIndex(fullaID,A,baseProp,K);
        Pz = accumarray(z,Pfull,[numel(aID) 1]);
        Nz = accumarray(z,1,[numel(aID) 1]);
        P = Pz(binID)./max(Nz(binID),1);
        
    elseif any(regexp(prop,'^expsmooth([A-Za-z]+)(_[0-9]+|$)'))
        %Exponentially smoothing from birth
        baseProp = regexp(prop,'^expsmooth([A-Za-z]+)(_[0-9]+|$)','tokens');
        alpha = baseProp{1}{2}; %smoothing parameter
        baseProp = baseProp{1}{1};
                
        if isempty(alpha)
            alpha = 0.5;
        else
            alpha = str2double(alpha(2:end))*0.01;
        end
        
        %Get the full cell cycle indices for each aID
        cell2aid = obj.cell2aid(A,K,false);
        
        %Find unique cellIDs corresponding to aID
        cellIDs = A(aID(:),obj.a('cellID'));
        
        %Get all relevant aID
        fullaID = cell2aid(cellIDs);
        
        %Truncate each to include the birth aID up to the aID of interest
        shortaID = cellfun(@(x,y) x(1:find(x==y)),fullaID,num2cell(aID(:)),'uniformoutput',0);
        
        %Retrieve the property of interest
        nID = cellfun(@numel,shortaID);
        Pc = obj.getPropAIndex(cell2mat(shortaID),A,baseProp,K);
        Pc = mat2cell(Pc,nID,1);
       
        %Perform exponential smoothing
        %First pad single aID cells so that at least two indices are
        %present
        doubleMe = nID==1;
        g = @(x) x([1; 1]);
        Pc(doubleMe) = cellfun(g,Pc(doubleMe),'uniformoutput',false);
        
        f = @(x) sum(alpha.*(1-alpha).^((numel(x)-2:-1:0)').*x(2:end)) + (1-alpha).^(numel(x)-1).*x(1);
        P = cellfun(f,Pc);
        
    elseif any(regexp(prop,'^dbexpsmooth([A-Za-z]*?)(optim|$)'))
        %Double exponentially smoothing from birth
        baseProp = regexp(prop,'^dbexpsmooth([A-Za-z]*?)(optim|$)','tokens');
        doOptim = baseProp{1}{2}; %smoothing parameter
        baseProp = baseProp{1}{1};
                
        %Get the full cell cycle indices for each aID
        cell2aid = obj.cell2aid(A,K,false);
        
        %Find unique cellIDs corresponding to aID
        cellIDs = A(aID(:),obj.a('cellID'));
        
        %Get all relevant aID
        fullaID = cell2aid(cellIDs);
        
        %Truncate each to include the birth aID up to the aID of interest
        shortaID = cellfun(@(x,y) x(1:find(x==y)),fullaID,num2cell(aID(:)),'uniformoutput',0);
        
        %Retrieve the property of interest
        nID = cellfun(@numel,shortaID);
        Pc = obj.getPropAIndex(cell2mat(shortaID),A,baseProp,K);
        Pc = mat2cell(Pc,nID,1);
       
        %Perform exponential smoothing
        %First pad single aID cells so that at least two indices are
        %present
        doubleMe = nID==1;
        g = @(x) x([1; 1]);
        Pc(doubleMe) = cellfun(g,Pc(doubleMe),'uniformoutput',false);
        
        if isempty(doOptim)
            ab = obj.getSmoothingParams(baseProp);
            P = cellfun(@(x) dbexpsmoothMP_lastVal(x,ab),Pc);
        else
            P = cellfun(@(x) dbexpsmoothMP_lastVal(x,dbexpsmoothMP_optim(x)),Pc);    
        end
        
    elseif any(regexp(prop,'^dbexpsmooth2([A-Za-z]*?)(optim|$)'))
        %Double exponentially smoothing from birth, smoothing performed in
        %reverse first to determine the initial values to use.
        baseProp = regexp(prop,'^dbexpsmooth2([A-Za-z]*?)(optim|$)','tokens');
        doOptim = baseProp{1}{2}; %smoothing parameter
        baseProp = baseProp{1}{1};
                
        %Get the full cell cycle indices for each aID
        cell2aid = obj.cell2aid(A,K,false);
        
        %Find unique cellIDs corresponding to aID
        cellIDs = unique(A(aID(:),obj.a('cellID')));
        
        %Bin each aID to its cell
        cellBins = [cellIDs; cellIDs(end)+1];
        [~,~,binID] = histcounts(A(aID(:),obj.a('cellID')),cellBins);
        
        %Get all relevant aID in ascending order of cellIDs
        fullaID = cell2aid(cellIDs);
        nID = cellfun(@numel,fullaID);
        Pc = obj.getPropAIndex(cell2mat(fullaID),A,baseProp,K);
        Pc = mat2cell(Pc,nID,1);
        
        %Perform double exponential smoothing on each cell cycle, with
        %reverse
        %Double single-point time series so that the smoothing function
        %works properly
        doubleMe = nID==1;
        g = @(x) x([1; 1]);
        Pc(doubleMe) = cellfun(g,Pc(doubleMe),'uniformoutput',false);
        
        if isempty(doOptim)
            ab = obj.getSmoothingParams(baseProp);
            Pc = cellfun(@(x) dbexpsmoothMP2(x,ab),Pc,'uniformoutput',0);
        else
            Pc = cellfun(@(x) dbexpsmoothMP2(x,dbexpsmoothMP2_optim(x)),Pc,'uniformoutput',0);    
        end
        
        %Restructure Pc to same index as aID
        Pc = Pc(binID);
        fullaID = fullaID(binID);
                
        %Select only the index of interest within each smoothed cell cycle
        P = cellfun(@(x,y,z) x(find(y==z,1)),Pc,num2cell(aID),fullaID);
        
             
    elseif any(regexp(prop,'frac(\w+)_(\w+)'))
        %One property divided by the other, the first property can also be
        %1 to get the inverse of a property
        baseProp = regexp(prop,'frac(\w+)_(\w+)','tokens');
        baseProp1 = baseProp{1}{1};
        baseProp2 = baseProp{1}{2};
        
        if str2double(baseProp1)~=1
            P1 = obj.getPropAIndex(aID,A,baseProp1,K);
        else
            P1 = 1;
        end
        
        P2 = obj.getPropAIndex(aID,A,baseProp2,K);
        
        P = P1./P2;
            
    elseif ~isempty(regexp(prop,'added([A-Za-z])+([01])+','tokens'))
        %Volume, length, etc added since birth. Convention: 'addedV0',
        %'addedV0','addedL0','addedL1', where 0 or 1 denotes addition since
        %the first cell cycle point of this cell (0) or the last cell cycle
        %point of its parent (1)
        %Map all cell lines to reduce time indexing cell cycle points
        if ~obj.aIDmap(1).built
            obj.doaIDmap;
        end
                
        baseProp = regexp(prop,'added([A-Za-z])+([01])+','tokens');
        baseProp = baseProp{1}{1};
        Pv = obj.getPropAIndex(aID,A,baseProp,K);
        Pv = Pv(:);
        useParent = str2double(prop(end));
        subaID = obj.getEndCyclePtaID(K,aID,A,'first',useParent);
        Pv0 = nan(numel(aID),1);
        Pv0(~isnan(subaID)) = obj.getPropAIndex(subaID(~isnan(subaID)),A,baseProp,K);
        
        if obj.isExtensive(baseProp)
            %If an extensive property: halve parent values
            P = Pv - Pv0*(1/(1+useParent)); 
        else
            P = Pv - Pv0;
        end
        
    elseif ~isempty(regexp(prop,'diff([A-Za-z])+([0-9]+)(_[01]|$)','tokens'))
        %Difference between a property at the given aID, and n time points
        %around (within the same cell or its lineage, if far away enough).
        %Note that n = odd, will use floor(n/2) points back and ceil(n/2)
        %forward.
        
        baseProp = regexp(prop,'diff([A-Za-z])+([0-9]+)(_[01]|$)','tokens');
        nPts = str2double(baseProp{1}{2});
        useParent = baseProp{1}{3};
        if isempty(useParent)
            useParent = false;
        else
            useParent = logical(str2double(useParent(end)));
        end
        baseProp = baseProp{1}{1};
        
        %Perform calculation twice is useParent to get the correct average at the fork
        %between different daughter lineages since it was originally
        %written to select one daughter at a time.
        P = zeros(numel(aID),useParent+1);
        for dd = 1:useParent+1
            %Get indices for all data points to calculate the difference
            D = obj.getdiffaID(K,aID,A,nPts,dd);
            
            %Retrieve factor for multiplication of e.g. extensive variables
            %or to remove other generations where requested
            [otherGenUsed,D] = obj.formatOtherGenUsed(D,useParent,obj.isExtensive(baseProp));
            
            daID = D.daID;
            ptRange = D.ptRange;
            
            %Get set of valid point ranges 
            ptRangeSet = unique(ptRange,'rows');
            ptRangeSet(ptRangeSet(:,1)==0,:) = [];
            
            %Get property of interest
            isValidIndex = ~isnan(daID);
            
            P0 = obj.getPropAIndex(daID(isValidIndex),A,baseProp,K);
            
            dP = nan(numel(aID),nPts+1);
            dP(isValidIndex) = P0;
            dP = dP.*otherGenUsed;
            
            if nPts==1
                %If only two time points are used (nPts = 1) just take the
                %difference
                P(:,dd) = diff(dP,1,2);
            else
                %Otherwise, perform a linear fit. Note that here, we assume
                %the time difference between points is approximately linear.
                %Subset of dP that has tracked values
                for j = 1:size(ptRangeSet,1)
                    useMe = ptRange==ones(size(ptRange,1),1)*ptRangeSet(j,:);
                    useMe = all(useMe,2);
                    
                    %Linear fit to get difference: assumes constant time step
                    xx = ptRangeSet(j,1):ptRangeSet(j,2);
                    xx = [xx(:) ones(numel(xx),1)];
                    yy = dP(useMe,ptRangeSet(j,1):ptRangeSet(j,2))';
                    pp = xx\yy;
                    
                    dPj = pp(1,:)*ptRangeSet(j,2) - pp(1,:)*ptRangeSet(j,1);
                    dPj = dPj./diff(ptRangeSet(j,:));
                    P(useMe,dd) = dPj;
                end
            end
        end
        
        P = mean(P,2);
        
        
    elseif ~isempty(regexp(prop,'^diff2(\w*?)_([1-9]+)(_[01])','tokens'))
        %Difference between a property at the given aID, and n time points
        %around (within the same cell or its lineage, if far away enough).
        %Note that n = odd, will use floor(n/2) points back and ceil(n/2)
        %forward.
        
        baseProp = regexp(prop,'diff2(\w*?)_([1-9]+)(_[01])','tokens');
        nPts = str2double(baseProp{1}{2});
        useParent = baseProp{1}{3};
        if isempty(useParent)
            useParent = false;
        else
            useParent = logical(str2double(useParent(end)));
        end
        baseProp = baseProp{1}{1};
        
        %Perform calculation twice to get the correct average at the fork
        %between different daughter lineages since it was originally
        %written to select only a random daughter
        P = zeros(numel(aID),2);
        for dd = 1:2
            %Get indices for all data points to calculate the difference
            D = obj.getdiffaID(K,aID,A,nPts,dd);
            
            %Retrieve factor for multiplication of e.g. extensive variables
            %or to remove other generations where requested
            [otherGenUsed,D] = obj.formatOtherGenUsed(D,useParent,obj.isExtensive(baseProp));
            
            daID = D.daID;
            ptRange = D.ptRange;
            
            %Get set of valid point ranges 
            ptRangeSet = unique(ptRange,'rows');
            ptRangeSet(ptRangeSet(:,1)==0,:) = [];
            
            %Get property of interest
            isValidIndex = ~isnan(daID);
            
            P0 = obj.getPropAIndex(daID(isValidIndex),A,baseProp,K);
            dP = nan(numel(aID),nPts+1);
            dP(isValidIndex) = P0;
            dP = dP.*otherGenUsed;
            
            if nPts==1
                %If only two time points are used (nPts = 1) just take the
                %difference
                P(:,dd) = diff(dP,1,2);
            else
                %Otherwise, perform a linear fit. Note that here, we assume
                %the time difference between points is approximately linear.
                %Subset of dP that has tracked values
                for j = 1:size(ptRangeSet,1)
                    useMe = ptRange==ones(size(ptRange,1),1)*ptRangeSet(j,:);
                    useMe = all(useMe,2);
                    
                    %Linear fit to get difference: assumes constant time step
                    xx = ptRangeSet(j,1):ptRangeSet(j,2);
                    xx = [xx(:) ones(numel(xx),1)];
                    yy = dP(useMe,ptRangeSet(j,1):ptRangeSet(j,2))';
                    pp = xx\yy;
                    
                    dPj = pp(1,:)*ptRangeSet(j,2) - pp(1,:)*ptRangeSet(j,1);
                    dPj = dPj./diff(ptRangeSet(j,:));
                    P(useMe,dd) = dPj;
                end
            end
        end
        
        P = mean(P,2);
        
    elseif ~isempty(regexp(prop,'^growthRate([A-Za-z])+([1-9]+)(_[01]|$)','tokens'))
        %Growth rate of a certain property using a particular number of
        %points        
        baseProp = regexp(prop,'^growthRate([A-Za-z])+([1-9]+)(_[01]|$)','tokens');
        nPts = str2double(baseProp{1}{2});
        useParent = baseProp{1}{3};
        baseProp = baseProp{1}{1};
                        
        dP0 = obj.getPropAIndex(aID,A,sprintf('diff%s%d%s',baseProp,nPts,useParent),K);
        dt = obj.getPropAIndex(aID,A,sprintf('difft%d%s',nPts,useParent),K);
        P0 = obj.getPropAIndex(aID,A,baseProp,K);
                
        P = dP0./dt./P0;
        
    elseif ~isempty(regexp(prop,'^growthRate2(\w*?)_([1-9])_([01])','tokens'))
        %Growth rate of a certain property using a particular number of
        %points
        baseProp = regexp(prop,'^growthRate2(\w*?)_([1-9])_([01])','tokens');
        nPts = str2double(baseProp{1}{2});
        useParent = baseProp{1}{3};
        baseProp = baseProp{1}{1};
        
        dP0 = obj.getPropAIndex(aID,A,sprintf('diff2%s_%d_%s',baseProp,nPts,useParent),K);
        dt = obj.getPropAIndex(aID,A,sprintf('difft%d_%s',nPts,useParent),K);
        P0 = obj.getPropAIndex(aID,A,baseProp,K);
        
        P = dP0./dt./P0;
        
    elseif ~isempty(regexp(prop,'^beta([A-Za-z]+)_([A-Za-z]+)([0-9]+)(_[01]|$)','tokens'))
        %Growth rate of a certain property, assumed to be proportional to the second property, e.g.
        %for beta in Harris-Theriot assumption
        baseProp = regexp(prop,'^beta([A-Za-z]+)_([A-Za-z]+)([0-9]+)(_[01]|$)','tokens');
        nPts = str2double(baseProp{1}{3});
        useParent = baseProp{1}{4};
        baseProp2 = baseProp{1}{2};
        baseProp1 = baseProp{1}{1};
                                
        dP0 = obj.getPropAIndex(aID,A,sprintf('diff%s%d%s',baseProp1,nPts,useParent),K);
        dt = obj.getPropAIndex(aID,A,sprintf('difft%d%s',nPts,useParent),K);
        P0 = obj.getPropAIndex(aID,A,baseProp2,K);
                
        P = dP0./dt./P0;
        
    elseif ~isempty(regexp(prop,'^beta2([a-zA-Z]\w*?)_([a-zA-Z]\w*?)_([0-9]+)_([01])','tokens'))
        %Growth rate of a certain property, assumed to be proportional to the second property, e.g.
        %for beta in Harris-Theriot assumption
        baseProp = regexp(prop,'^beta2([a-zA-Z]\w*?)_([a-zA-Z]\w*?)_([0-9]+)_([01])','tokens');
        nPts = str2double(baseProp{1}{3});
        useParent = baseProp{1}{4};
        baseProp2 = baseProp{1}{2};
        baseProp1 = baseProp{1}{1};
                        
        dP0 = obj.getPropAIndex(aID,A,sprintf('diff2%s_%d_%s',baseProp1,nPts,useParent),K);
        dt = obj.getPropAIndex(aID,A,sprintf('difft%d_%s',nPts,useParent),K);
        P0 = obj.getPropAIndex(aID,A,baseProp2,K);
                      
        P = dP0./dt./P0;
                
    elseif ~isempty(regexp(prop,'^loggrowthRate([A-Za-z])+([0-9]+)(_[01]|$)','tokens'))
        %Growth rate of a certain property using a particular number of
        %points calculated by a linear fitting of the log value, with
        %option to exclude or include parent or daughter cell cycle points
        baseProp = regexp(prop,'^loggrowthRate([A-Za-z])+([0-9]+)(_[01]|$)','tokens');
        nPts = str2double(baseProp{1}{2});
        useParent = baseProp{1}{3};
        if isempty(useParent)
            useParent = false;
        else
            useParent = logical(str2double(useParent(end)));
        end
        baseProp = baseProp{1}{1};
        
        %Perform calculation twice to get the correct average at the fork
        %between different daughter lineages 
        P = nan(numel(aID),2);
        for dd = 1:2
            %Get indices for all data points to calculate the difference
            D = obj.getdiffaID(K,aID,A,nPts,dd);
            
            %Factor for multiplication to account for cross generation
            %calculations or remove parents/daughters if requested
            [otherGenUsed,D] = obj.formatOtherGenUsed(D,useParent,obj.isExtensive(baseProp));
            
            daID = D.daID;
            ptRange = D.ptRange;
            
            %Get set of valid point ranges
            ptRangeSet = unique(ptRange,'rows');
            ptRangeSet(ptRangeSet(:,1)==0,:) = [];
            
            %Get property of interest
            isValidIndex = ~isnan(daID);
            
            P0 = obj.getPropAIndex(daID(isValidIndex),A,baseProp,K);
            dP = nan(numel(aID),nPts+1);
            dP(isValidIndex) = P0;
            dP = dP.*otherGenUsed;
            
            t0 = obj.getPropAIndex(daID(isValidIndex),A,'t',K);
            dt = nan(numel(aID),nPts+1);
            dt(isValidIndex) = t0;
            dt = diff(dt,1,2);
            dt = nanmean(dt(:));
            dt = (0:nPts)*dt;
            
            for j = 1:size(ptRangeSet,1)
                if ptRangeSet(j,1)==ptRangeSet(j,2)
                    continue
                end
                useMe = ptRange==ones(size(ptRange,1),1)*ptRangeSet(j,:);
                useMe = all(useMe,2);

                %Linear fit to get difference: assumes constant time step
                xx = dt(ptRangeSet(j,1):ptRangeSet(j,2));
                xx = [xx(:) ones(numel(xx),1)];
                yy = log(dP(useMe,ptRangeSet(j,1):ptRangeSet(j,2))');
                pp = xx\yy;
                P(useMe,dd) = pp(1,:);
            end
        end
        
        P = mean(P,2);
        
    elseif ~isempty(regexp(prop,'d([A-Za-z])+byd([A-za-z])+([1-9])+(_[01]|$)','tokens'))
        %First derivatives dX by dY using n points
        baseProp = regexp(prop,'d([A-Za-z])+byd([A-za-z])+([1-9])+(_[01]|$)','tokens');
        baseProp1 = baseProp{1}{1};
        baseProp2 = baseProp{1}{2};
        nPts = str2double(baseProp{1}{3});
        useParent = baseProp{1}{4};
        
        dP1 = obj.getPropAIndex(aID,A,sprintf('diff%s%d%s',baseProp1,nPts,useParent),K);
        dP2 = obj.getPropAIndex(aID,A,sprintf('diff%s%d%s',baseProp2,nPts,useParent),K);
        
        P = dP1./dP2;
        
    elseif any(regexp(prop,'^smooth\w+_[0-9]+'))
        %Smoothed with moving average along cell cycle using a window of N
        %points, with ceil((N-1)/2) points forward and floor((N-1)/2)
        %points back. No daughter/mother cell points are used
        baseProp = regexp(prop,'smooth(\w+)_([0-9]+)','tokens');
        nPts = str2double(baseProp{1}{2});
        baseProp = baseProp{1}{1};
        
        %Get the data points for smoothing
        D = obj.getdiffaID(K,aID,A,nPts,1);
        
        %Format otherGenUsed to exclude out of cycle points 
        [otherGenUsed,D] = obj.formatOtherGenUsed(D,false,false);
        daID = D.daID;
        
        %Get property of interest
        isValidIndex = ~isnan(daID);
        
        P0 = obj.getPropAIndex(daID(isValidIndex),A,baseProp,K);
        dP = nan(numel(aID),nPts+1);
        dP(isValidIndex) = P0;
        dP = dP.*otherGenUsed;
        
        P = nanmean(dP,2);
        
    elseif any(regexp(prop,'^expfitA[a-zA-Z]+[0-9]+_[01]'))
        %Exponential fit using a sliding window of nPts, with the option to
        %use or ignore points from other generations
        baseProp = regexp(prop,'^expfitA([a-zA-Z]+)([0-9]+)_([01])','tokens');
        nPts = str2double(baseProp{1}{2});
        useParent = str2double(baseProp{1}{3});
        baseProp = baseProp{1}{1};
        
        %Perform calculation twice to get the correct average at the fork
        %between different daughter lineages
        P = nan(numel(aID),useParent+1);
        for dd = 1:useParent+1
            %Get indices for all data points to calculate the difference
            D = obj.getdiffaID(K,aID,A,nPts,dd);
            
            %Factor for multiplication to account for cross generation
            %calculations or remove parents/daughters if requested
            [otherGenUsed,D] = obj.formatOtherGenUsed(D,useParent,obj.isExtensive(baseProp));
            
            daID = D.daID;
            daID(isnan(daID)) = 0;
            
            ptRange = D.ptRange;
            
            %Get set of valid point ranges
            ptRangeSet = unique(ptRange,'rows');
            ptRangeSet(ptRangeSet(:,1)==0,:) = [];
            
            %Get property of interest            
            dP = obj.getPropAIndex(daID,A,baseProp,K);
            dP = dP.*otherGenUsed;
            
            for j = 1:size(ptRangeSet,1)
                useMe = ptRange==ones(size(ptRange,1),1)*ptRangeSet(j,:);
                useMe = all(useMe,2);
                
                if ptRangeSet(j,1)==ptRangeSet(j,2)
                    P(useMe,dd) = dP(useMe,ptRangeSet(j,1));
                else
                    %Linear fit to get difference: assumes constant time step
                    xx = ptRangeSet(j,1):ptRangeSet(j,2);
                    xx = [xx(:) ones(numel(xx),1)];
                    yy = log(dP(useMe,ptRangeSet(j,1):ptRangeSet(j,2))');
                    pp = xx\yy;

                    P(useMe,dd) = exp(pp(1,:).*(floor(nPts/2)+1) + pp(2,:));
                end
            end
        end
        P = mean(P,2);
        
    elseif ~isempty(regexp(prop,'subtract(\w+)from(\w+)','tokens'))
        baseProp = regexp(prop,'subtract(\w+)from(\w+)','tokens');
        P1 = obj.getPropAIndex(aID(:),A,baseProp{1}{1},K);
        P2 = obj.getPropAIndex(aID(:),A,baseProp{1}{2},K);
        
        P = P2-P1;
        
    elseif regexp(prop,'^SAareadbexpsm')
        %Surface area calculated from double exponentially smoothed length
        %and cross-sectional area
        L = obj.getPropAIndex(aID(:),A,'dbexpsmooth2L',K);
        w = obj.getPropAIndex(aID(:),A,'wareadbexpsm',K);
        area = A(aID,obj.a('area'));
        
        P = obj.surfacearea(L,w,area,'hemisph');
        
    elseif regexp(prop,'^SA')
        %Surface area
        L = obj.getPropAIndex(aID(:),A,'L',K);
        wRaw = obj.getPropAIndex(aID(:),A,'w',K);
        w = wRaw;
        area = A(aID,obj.a('area'));
        
        P = obj.surfacearea(L,w,area,prop(3:end));
              
    elseif any(regexp(prop,'tUntilDiv([0-9]+|$)(_[01]|$)'))
        %Time until Nth division, where N is the generation number. N>=0,
        %and NaNs returned if the Nth generation has not been tracked. Note
        %that if N > 0, the returned output P will be a matrix where rows
        %correspond to the provided aID with 2^N columns for each possible
        %lineage. Last option to use the current data filter to exclude
        %some results
        
        %Determine generation number
        nGen = regexp(prop,'tUntilDiv([0-9]+|$)(_[01]|$)','tokens');
        useFilter = nGen{1}{2};
        nGen = nGen{1}{1};
        if isempty(nGen)
            nGen = 0;
            useFilter = false;
        else
            nGen = str2double(nGen);
            if isempty(useFilter)
                useFilter = true;
            else
                useFilter = logical(str2double(useFilter(end)));
            end
        end
        
        if useFilter
            [fov,L] = obj.fovLineK(K);
        end
        
        %Get index of zeroth division
        firstLast = 'last';
        
        if nGen==0
            aID1 = obj.getEndCyclePtaID(K,aID,A,firstLast,false);
            P = A(aID1,obj.a('t')) - A(aID,obj.a('t'));
        else
            aID1 = zeros(numel(aID0),2^nGen);
            aID1k = aID(:);
            for k = 1:nGen
                aID1k1 = obj.getEndCyclePtaID(K,aID1k,A,'last',1,1);
                aID1k2 = obj.getEndCyclePtaID(K,aID1k,A,'last',1,2);
                aID1k = [aID1k1 aID1k2]';
                aID1k = aID1k(:);
                
                if useFilter
                    bID = obj.getBIndex(fov,L,obj.getPropAIndex(aID1k,A,'cellID',K));
                    passFilter = obj.getPropBIndex(bID,'passFilter');
                    aID1k(not(passFilter>0)) = 0;
                end
            end
            aID1k = obj.getEndCyclePtaID(K,aID1k,A,'last',0,0);
            aID1k = reshape(aID1k,2^nGen,numel(aID(:)))';
            
            %Get only unique daughter points
            isUnique = aID1k(:,1)~=aID1k(:,2);
            aID1k(~isUnique,2) = 0;
            
            %Return to matrix of appropriate size
            aID1(isfinite(aID0(:)) & (aID0(:)~=0),:) = aID1k;
            aID1(isnan(aID1)) = 0;
            
            %Get time until division
            P = repmat(obj.getPropAIndex(aID0(:),A,'t',K),1,2^nGen);
            P = obj.getPropAIndex(aID1,A,'t',K) - P;
            return
        end
        
    elseif strcmp(prop,'cellPhase')
        %Phase during cell cycle
        aID0 = obj.getEndCyclePtaID(K,aID,A,'first',0); %birth index
        aID1 = obj.getEndCyclePtaID(K,aID,A,'last',0);  %division index
         
        P = 1 - (A(aID1,obj.a('t')) - A(aID,obj.a('t')))./(A(aID1,obj.a('t')) - A(aID0,obj.a('t')));
                
    elseif strcmp(prop,'tTrans')
        %Time translated by shift of interest
        P = A(aID,obj.a('t')) - obj.currShift;        
    else
        %Retrieve property from a column of the line matrices
        P = A(aID,obj.a(prop));
    end
    
    %Convert growth rates to doublings per hour
    if any(strfind(prop,'dblH'))
        P = P.*60./log(2);
    elseif any(strfind(prop,'perH'))
        %Convert to per hour (default is per minute)
        P = P.*60;
    elseif any(strfind(prop,'inH'))
        %Convert time to hours
        P = P./60;
    end
    
    %Return in same shape as aID
    PaID(aID0~=0 & isfinite(aID0)) = P;
    P = PaID;
    
    end
    
    
    
    
    
    function P = getPropBIndex(obj,bID,prop)
    %P = obj.getPropBIndex(bID,prop)
    %Retrieves the property prop from matrix B according to indices in B.
    %If the requested indices from B, bID is not a logical, P is returned 
    %in the same shape as bID, with NaNs where no cell was found; otherwise
    %B is returned as a vector with the appropriate number of elements and
    %NaNs where no property was found.
    
    if ~islogical(bID)
        P0 = nan(size(bID));
        bID0 = bID;
        bID = bID(bID~=0 & isfinite(bID));
    else
        bID0 = bID;
        bID = find(bID);
    end
        
    if strcmp(prop,'t1')
        %Time of division
        %Determine the lines and A indices that need to be read
        [K,aID] = obj.bID2aID(bID);
        
        P = nan(numel(bID),1);
        for k = unique(K')
            A = obj.loadLineK(k);
            aIDk = aID(K==k);
            nID = cellfun(@numel,aIDk);
            aIDk = cell2mat(aIDk);
            Pk = obj.getPropAIndex(aIDk,A,'t',k);
            Pk = mat2cell(Pk,nID,1);
            P(K==k) = cellfun(@(x) x(end),Pk);
        end
        
    elseif regexp(prop,'^t[01]Trans$')
        %Time translated by the shift time
        timePt = regexp(prop,'^t([01])Trans$','tokens');
        timePt = timePt{1}{1};
        
        P = obj.getPropBIndex(bID(:),['t' timePt]) - obj.currShift;
        
    elseif strcmp(prop,'frame1')
        %Frame at which division occurrred
        %Determine the lines and A indices that need to be read
        [K,aID] = obj.bID2aID(bID);
        nID = cellfun(@numel,aID);
        
        P = nan(numel(bID),1);
        for k = unique(K')
            A = obj.loadLineK(k);
            aIDk = aID(K==k);
            nIDk = nID(K==k);
            aIDk = vertcat(aIDk{:});
            Pk = obj.getPropAIndex(aIDk,A,'frame',k);
            P(K==k) = Pk(cumsum(nIDk)); 
        end
                
    elseif strcmp(prop(1),'w')
        %Derive from area and length
        area = obj.getPropBIndex(bID(:),['area' prop(2)]);
        L = obj.getPropBIndex(bID(:),['L' prop(2)]);
                
        P = 2*(L-sqrt(L.^2-(4-pi)*area))./(4-pi);
        
    elseif any(regexp(prop,'^smth[a-zA-z]+[01]_[0-9]+'))
        %Cell-cycle smoothed property at either the beginning or end of the
        %cell cycle
        baseProp = regexp(prop,'^smth([a-zA-z]+)([01])_([0-9]+)','tokens');
        firstLast = not(str2double(baseProp{1}{2}));
        baseProp = sprintf('smooth%s_%s',baseProp{1}{1},baseProp{1}{3});
        
        [K,aID] = obj.bID2aID(bID(:));
        nID = cellfun(@numel,aID);
        
        P = nan(numel(bID),1);
        for k = unique(K');
            A = obj.loadLineK(k);
            aIDk = aID(K==k);
            aIDk = vertcat(aIDk{:});
            nIDk = nID(K==k);
            
            idk = cumsum(nIDk) + firstLast;
            idk = idk(1:end-firstLast);
            idk = [firstLast; idk]; %#ok<AGROW>
            idk = idk(idk~=0);
            
            aIDk = aIDk(idk);
            Pk = obj.getPropAIndex(aIDk,A,baseProp,k);
            P(K==k) = Pk;
        end
        
    elseif any(regexp(prop,'^expfitB[a-zA-Z]+[0123]'))
        %Exponential fit for either the beginning (0) or end (1) of the
        %cell cycle. Use option 2 to retrieve the exponent in units
        %min^-1, option 3 to retrieve the prefactor
        baseProp = regexp(prop,'^expfitB([a-zA-Z]+)([0123])','tokens');
        firstLast = str2double(baseProp{1}{2});
        baseProp = baseProp{1}{1};
        
        [K,aID] = obj.bID2aID(bID(:));
        nID = cellfun(@numel,aID);
        
        %If growth rate requested, convert to min^{-1}
        do_convert = firstLast==2;
        if do_convert
            dt = cell(numel(bID),1);
        end
        
        %Retrieve cell cycle properties for each cell
        Pc = cell(numel(bID),1);
        for k = unique(K')
            A = obj.loadLineK(k);
            aIDk = aID(K==k);
            nIDk = nID(K==k);
            Pk = obj.getPropAIndex(vertcat(aIDk{:}),A,baseProp,k);
            Pk = mat2cell(Pk,nIDk,1);
            Pc(K==k) = Pk;
            
            %Retrieve dt to get frames per min for conversion
            if do_convert
                dtk = obj.getPropAIndex(vertcat(aIDk{:}),A,'difft1',k);
                dtk = mat2cell(dtk,nIDk,1);
                dt(K==k) = dtk;
            end
        end
        
        %Perform exponential fitting by linear regression on log
        P = nan(numel(bID),1);
        for n = unique(nID')
            xx = [(1:n)' ones(n,1)];
            yy = Pc(nID==n);
            
            if n~=1
                yy = log([yy{:}]);
                pp = xx\yy;
                Pn = [exp([pp(1,:).*0.5; pp(1,:).*(n+0.5)] + pp([2 2],:)); pp(1,:); exp(pp(2,:))];
            else
                Pn = [[yy{:}]; [yy{:}]; nan(2,numel(yy))];
            end
                                        
            P(nID==n) = Pn(firstLast+1,:);
        end
        
        %Perform conversion to min{-1} if necessary
        if do_convert
            dt_med = nanmedian(cell2mat(dt));
            P = P./dt_med;
        end
                
        %Perform exponential fitting by linear regression on log
        P = nan(numel(bID),1);
        for n = unique(nID')
            xx = [(1:n)' ones(n,1)];
            yy = Pc(nID==n);
            
            if n~=1
                yy = log([yy{:}]);
                pp = xx\yy;
                Pn = [exp([pp(1,:).*0.5; pp(1,:).*(n+0.5)] + pp([2 2],:)); pp(1,:); exp(pp(2,:))];
            else
                Pn = [[yy{:}]; [yy{:}]; nan(2,numel(yy))];
            end
                                        
            P(nID==n) = Pn(firstLast+1,:);
        end
        
        %Perform conversion to min{-1} if necessary
        if do_convert
            dt_med = nanmedian(cell2mat(dt));
            P = P./dt_med;
        end
        
                
    elseif strcmp(prop(1),'L')
        %Add pixel offset to length
        P = obj.flatData.B(bID,obj.b(prop)) + 0.1067;
        
    elseif regexp(prop,'^mean\w')
        %The mean cell cycle value of particular property
        baseProp = regexp(prop,'mean(\w+)','tokens');
        baseProp = baseProp{1}{1};
        
        %Determine the lines and A indices that need to be read
        [K,aID] = obj.bID2aID(bID);
                
        P = nan(numel(bID),1);
        for k = unique(K')
            A = obj.loadLineK(k);
            aIDk = aID(K==k);
            nID = cellfun(@numel,aIDk);
            aIDk = cell2mat(aIDk);
            Pk = obj.getPropAIndex(aIDk,A,baseProp,k);
            Pk = mat2cell(Pk,nID,1);
            P(K==k) = cellfun(@nanmean,Pk);
        end
        
    elseif strcmp(prop(1),'V')
        %VOLUME
        %Format as Vxyz01, where xyz is the method for volume calculation,
        %and 01 is 0 for the first and 1 for the last cell cycle points
        baseProp = regexp(prop,'V([a-z])+([01])','tokens');
        if ~isempty(baseProp)
            firstLast = baseProp{1}{2};
            mthd = baseProp{1}{1};
        elseif regexp(prop,'V[01]')
            %If no method specified, defaults to hemispherical caps
            firstLast = prop(2);
            mthd = 'hemisph';
        end
        
        L = obj.getPropBIndex(bID,['L' firstLast]);
        w = obj.getPropBIndex(bID,['w' firstLast]);
        area = obj.getPropBIndex(bID,['area' firstLast]);
        
        P = obj.volume(L,w,area,mthd);
                
        
    elseif regexp(prop,'^SA') %any(strfind(prop,'SA')) && isempty(strfind(prop,'growthRate'));
        %SURFACE AREA
        %Format as SAxyz01, where xyz is the method for calculation,
        %and 01 is 0 for the first and 1 for the last cell cycle points
        baseProp = regexp(prop,'^SA([a-z])+([01])','tokens');
        if ~isempty(baseProp)
            firstLast = baseProp{1}{2};
            mthd = baseProp{1}{1};
        elseif regexp(prop,'SA[01]')
            %If no method specified, defaults to hemispherical caps
            firstLast = prop(end);
            mthd = 'hemisph';
        end
        
        L = obj.getPropBIndex(bID,['L' firstLast]);
        w = obj.getPropBIndex(bID,['w' firstLast]);
        area = obj.getPropBIndex(bID,['area' firstLast]);
        
        P = obj.surfacearea(L,w,area,mthd);
        
        
    elseif regexp(prop,'^C([a-z]+[01]|[01])')
        %CONCENTRATION
        %Format as Cxyz01, where xyz is the method used to calculate
        %volume, and 01 is 0 for the first point of the cell cycle and 1
        %for the last
        baseProp = regexp(prop,'^C([a-z])+([01])','tokens');
        if ~isempty(baseProp)
            firstLast = baseProp{1}{2};
            mthd = baseProp{1}{1};
        elseif regexp(prop,'C[01]')
            %If no method specified, defaults to hemispherical caps
            firstLast = prop(2);
            mthd = 'hemisph';
        else
            warning('Invalid format for concentration. Empty matrix returned. \n')
            P = [];
            return
        end
        
        L = obj.getPropBIndex(bID,['L' firstLast]);
        w = obj.getPropBIndex(bID,['w' firstLast]);
        area = obj.getPropBIndex(bID,['area' firstLast]);
        
        V = obj.volume(L,w,area,mthd);
            
        P = obj.flatData.B(bID,obj.b(['Itot' firstLast]))./V(:);
        
    elseif regexp(prop,'delta[a-zA-z]+')
        %DIFFERENCE BETWEEN FIRST AND LAST CELL CYCLE POINTS 
        baseProp = regexp(prop,'delta([a-zA-z]+)','tokens');
        baseProp = baseProp{1}{1};
        
        P = obj.getPropBIndex(bID,[baseProp '1']) - obj.getPropBIndex(bID,[baseProp '0']);
        
    elseif regexp(prop,'(alpha[LV])')
        %Growth rate (determined from the entire cell cycle)
        baseProp = regexp(prop,'(alpha[LV])','tokens');
        baseProp = baseProp{1}{1};
        
        P = obj.flatData.B(bID,obj.b(baseProp));
                
    elseif strcmp(prop,'passFilter')
        %Check if the cell of interest passes the current filter
        P = obj.filterIDFlat(bID);
        
    elseif isKey(obj.b,prop)
        %Anything other property that is contained in the flat file
        P = obj.flatData.B(bID,obj.b(prop));    
    else
        %No property retrieved
        fprintf('WARNING: %s is not a valid property to request from B matrix data, empty matrix returned.\n',prop)
        P = [];
        return
    end
    
%     %Translate time if needed
%     if any(strfind(prop,'Trans'))
%         P = P - obj.currShift;
%     end
    
    %Convert times if needed
    if any(strfind(prop,'dblH'))
        P = P.*60./log(2);
    elseif any(strfind(prop,'perH'))
        P = P.*60;
    elseif any(strfind(prop,'inH'))
        P = P./60;
    end
        
    if ~islogical(bID0)
        P0(bID0~=0 & isfinite(bID0)) = P;
    else
        P0 = P;
    end
    
    P = P0;
    
    end
    
    
    
    function [propA,propB] = propParse(obj,prop)
    %[propA,propB] = obj.propParse(prop)
    %Separates requested properties into two arrays, one for properties
    %that should be retrieved using A (i.e. instantaneous properties) and 
    %another for B (i.e. flat properties) using try/catch filter.
    
    if ~iscell(prop); prop = {prop}; end
    
    propA = {};
    ac = 1;
    propB = {};
    bc = 1;
    rejected = {};
    rc = 1;
    
    %Load dummy A index data
    A = obj.loadLineK(1);
    aID = (1:5)';
    
    for p = 1:numel(prop)
        try
            p0 = obj.getPropBIndex(1,prop{p});
            if isempty(p0)
                error %#ok<LTARG>
            end
            propB{bc} = prop{p};
            bc = bc + 1;
        catch
            try 
               p0 = obj.getPropAIndex(aID,A,prop{p},1);
               propA{ac} = prop{p};
               ac = ac + 1;
            catch
               rejected{rc} = prop{p};
               rc = rc + 1;
            end
        end
    end
    
    if numel(rejected)>0
        fprintf('Unable to retrieve property: %s \n',rejected{:})
    end
    
    end
    
    
    
    function D = getdiffaID(obj,K,aID,A,nPts,daughter2use)
    %D = obj.getdiffaID(obj,K,aID,A,nPts,daughter2use). Retrieves the indices up to and 
    %including aID +/- nPts/2, as determined by cell lineage. Used by 'diff' 
    %option of getPropAIndex e.g. for the calculation of derivatives. The 
    %output is stored in the object (diffaID) to save time in future calculations.
    %D is a struct with fields:
    %daID           numel(aID) by nPts+1 (since it includes each aID, the point of interest)
    %               indices in A corresponding to +/- nPts/2 from each aID
    %
    %parentUsed     logical indicating if the index refers to a parental cell
    %               note that the last column simply corresponds to aID,
    %               therefore it will always consist of zeros
    %
    %daughterUsed   logical indicating if the index refers to a daughter cell
    %
    %ptRange        numel(aID) by 2 indicating the range of indices in daID to
    %               use since NaNs are in daID where no cell was tracked
    %               (saves time that would be spent searching for valid indices later)
    %

    D = struct;

    try
        %Check if it has been calculated and stored in the object already
        daID = obj.diffaID(K,nPts).daID{daughter2use}(aID,:);
        parentUsed = obj.diffaID(K,nPts).parentUsed{daughter2use}(aID,:);
        daughterUsed = obj.diffaID(K,nPts).daughterUsed{daughter2use}(aID,:);
        ptRange = obj.diffaID(K,nPts).ptRange{daughter2use}(aID,:);
        
        if any(sum(daID,2)==0) || any(all(isnan(daID),2))
            error('Some invalid daID retrieved. Attempting to recalculate.')
        end
    catch    
        %Get indices for all data points to calculate the difference
        daIDrev = zeros(numel(aID),floor(nPts/2)+1);
        daIDrev(:,end) = aID;
        parentUsed = false(numel(aID),floor(nPts/2)+1);
        parentFirstNaN = ones(numel(aID),1);
        for j = floor(nPts/2):-1:1
            [aIDj,parentUsedj] = obj.getNextCyclePtaID(K,daIDrev(:,j+1),A,'prev');
            daIDrev(:,j) = aIDj;
            parentUsed(:,j) = parentUsedj;
            parentFirstNaN(isnan(aIDj)) = max(parentFirstNaN(isnan(aIDj)),j+1);
        end

        daIDfor = zeros(numel(aID),ceil(nPts/2));
        daughterUsed = false(numel(aID),ceil(nPts/2));
        daughterFirstNaN = ones(numel(aID),1)*(nPts+1);
        for j = 1:ceil(nPts/2)
            if j==1
                [aIDj,daughterUsedj] = obj.getNextCyclePtaID(K,aID,A,'next',daughter2use);
            else
                [aIDj,daughterUsedj] = obj.getNextCyclePtaID(K,daIDfor(:,j-1),A,'next',daughter2use);
            end
            
            daIDfor(:,j) = aIDj;
            daughterUsed(:,j) = daughterUsedj;
            daughterFirstNaN(isnan(aIDj)) = min(daughterFirstNaN(isnan(aIDj)),j+floor(nPts/2));
        end

        daID = [daIDrev daIDfor];
        ptRange = [parentFirstNaN daughterFirstNaN];    
    end

    D.daID = daID;
    D.parentUsed = parentUsed;
    D.daughterUsed = daughterUsed;
    D.ptRange = ptRange;

    %Store in object
    obj.diffaID(K,nPts).daID{daughter2use}(aID,:) = daID;
    obj.diffaID(K,nPts).parentUsed{daughter2use}(aID,:) = parentUsed;
    obj.diffaID(K,nPts).daughterUsed{daughter2use}(aID,:) = daughterUsed;
    obj.diffaID(K,nPts).ptRange{daughter2use}(aID,:) = ptRange;

    end
    
    
    
    
    
    
    function [subaID,nextGenUsed] = getNextCyclePtaID(obj,K,aID,A,prevNext,daughter2use)
    %[subaID,nextGenUsed] = getNextCyclePtaID(obj,K,aID,A,prevNext,daughter2use)
    %Given a set of indices in A, aID, retrieves the index in A that
    %corresponds to either the previous or the next cell cycle point.
    %subaID is the same size as aID, with NaNs if the relevant cell cycle
    %point cannot be found. In the case of next cell cycle point requests,
    %a random daughter is used if the respective point in aID is the
    %division point; for previous points, the parent is used. nextGenUsed 
    %is a logical indicating where this is the case.
    %prevNext       set to 'prev' or 'next'
    
%     if obj.aIDmap(1).built
%         subaID = obj.aIDmap(K).(prevNext)(aID);
%         nextGenUsed = obj.aIDmap(K).([prevNext 'genUsed'])(aID);
%         
%         subaID = reshape(subaID,size(aID,1),size(aID,2));
%         nextGenUsed = reshape(nextGenUsed,size(aID,1),size(aID,2));
%         return
%     end
    
    aID = aID(:);
    
    subaID = nan(numel(aID),1);
    nextGenUsed = false(numel(aID),1);
    isValidID = (~isnan(aID)) & (aID~=0); %make subset if aID contains any NaNs
    
    if sum(isValidID)==0
        return
    end
    
    aID0 = aID(isValidID);
    subaID0 = nan(numel(aID0),1);
    nextGenUsed0 = false(numel(aID0),1);
    
    switch prevNext
        case 'prev'
            %Since data matrix A is written chronologically within a given cell,
            %can first see if aID - 1  is the correct index.
            aID_ = aID0 - 1;
            aID_nonzero = aID_;
            aID_nonzero(aID_<=0) = aID0(aID_<=0);
            
            sameCell = A(aID0,obj.a('cellID'))==A(aID_nonzero,obj.a('cellID'));
            sameCell = sameCell & (aID_>0);
            
            %For the rest, look for parents
            pt2use = 'first';
            
        case 'next'          
            %Since data matrix A is written chronologically within a given cell,
            %check aID + 1
            aID_ = aID0 + 1;
            aID_short = aID_;
            aID_short(aID_>size(A,1)) = aID0(aID_>size(A,1));
                        
            sameCell = A(aID0,obj.a('cellID'))==A(aID_short,obj.a('cellID'));
            sameCell = sameCell & (aID_<=size(A,1));
            
            %If next point is beyond a cell of interest, look for daughters
            pt2use = 'last';           
    end
    
    %If aID +/- 1 is within the same cell, store these
    subaID0(sameCell) = aID_(sameCell);
    
    %For the rest, look for parents/daughters
    if nargin < 6
        daughter2use = 0;
    end
    
    if any(~sameCell)
        otherGenaID = obj.getEndCyclePtaID(K,aID0(~sameCell),A,pt2use,1,daughter2use);
        subaID0(~sameCell) = otherGenaID;
        nextGenUsed0(~sameCell) = isfinite(otherGenaID);
    end
    
    %Store back into vector the same size as the original aID
    subaID(isValidID) = subaID0;
    nextGenUsed(isValidID) = nextGenUsed0;
    
    end
    
    
    function subaID = getEndCyclePtaID(obj,K,aID,A,firstLast,useNextGen,daughter2use)
    %subaID = getEndCyclePtaID(obj,K,aID,A,firstLast,useNextGen,daughter2use)
    %Given a set of indices in A, aID, retrieves in the index in A that
    %corresponds to either the first or last cell cycle point for the cell
    %respective to the index at aID. subaID is the same size as aID, with NaNs
    %if parent or daughter indices are requested but not found.
    %firstLast      set to 'first' or 'last' for respective cell cycle points
    %useNextGen     set to 1 (true) to use the last cycle point of the parent
    %                   in the case of 'first', or the first cell cycle point
    %                   of a random daughter in the case of 'last'. If set to 0
    %                   (false), only the present generation is used.
    
    subaID_ = zeros(size(aID));
    aID_ = aID;
    aID = aID(aID~=0 & isfinite(aID));
    
    if isempty(aID)
        subaID = subaID_;
        return
    end
    
    cellIDs = A(aID,[obj.a('cellID') obj.b('parentID')]);
    cellIDsUnique = unique(cellIDs,'rows');
    
    aID0 = zeros(size(cellIDsUnique,1),1);

    switch firstLast
        case 'first'
            if ~useNextGen
                %Get the first cell cycle point
                [aIDfull,cellBin] = ismember(A(:,obj.a('cellID')),cellIDsUnique(:,1));
                aIDfull = find(aIDfull);
                frame1 = logical([1; diff(A(aIDfull,obj.a('cellID')))]);
                frame1 = aIDfull(frame1);
                cellBin = cellBin(frame1);

                aID0(cellBin) = frame1;
            else
                %Get the last cell cycle point of each parent
                parentIDs = unique(cellIDsUnique(:,2));
                parentIDs(parentIDs==0) = [];
                [aIDfull,cellBin] = ismember(A(:,obj.a('cellID')),parentIDs);
                
                if sum(aIDfull)~=0
                    aIDfull = find(aIDfull);

                    frameEnd = logical([diff(A(aIDfull,obj.a('cellID'))); 1]);
                    frameEnd = aIDfull(frameEnd);
                    cellBin = cellBin(frameEnd);

                    map2child = zeros(max(parentIDs),1);
                    map2child(parentIDs(cellBin)) = frameEnd;

                    hasParent = cellIDsUnique(:,2)~=0;
                    aID0(hasParent) = map2child(cellIDsUnique(hasParent,2));
                end
            end

        case 'last'
            if ~useNextGen
                %Get the last cell cycle point
                [aIDfull,cellBin] = ismember(A(:,obj.a('cellID')),cellIDsUnique(:,1));
                aIDfull = find(aIDfull);
                frameEnd = logical([diff(A(aIDfull,obj.a('cellID'))); 1]);
                frameEnd = aIDfull(frameEnd);
                cellBin = cellBin(frameEnd);

                aID0(cellBin) = frameEnd;
            else
                %Get the first cell cycle point of each daughter
                [aIDfull,cellBin] = ismember(A(:,obj.a('parentID')),cellIDsUnique(:,1));
                
                if sum(aIDfull)~=0                
                    aIDfull = find(aIDfull);
                    frame1 = logical([1; diff(A(aIDfull,obj.a('cellID')))]);

                    frame1 = aIDfull(frame1);
                    cellBin = cellBin(frame1);

                    %Select the first, second, or a random daughter
                    if nargin < 7; daughter2use = 0; end

                    if daughter2use~=0
                        %Sort frame1 so that an ordered daughter is chosen
                        sortModes = {'descend','ascend'};

                        [frame1,sortID] = sort(frame1,1,sortModes{daughter2use});
                        cellBin = cellBin(sortID);
                    else
                        %Permute frame1 so that a random daughter can be chosen
                        rp = randperm(numel(frame1));
                        frame1 = frame1(rp);
                        cellBin = cellBin(rp);
                    end

                    aID0(cellBin) = frame1;
                end
            end
    end

    %Indices to build map back to aID
    edges = cellIDsUnique(:,1);
    edges = [edges; edges(end) + 1];
    [~,~,binID] = histcounts(cellIDs(:,1),edges);

    subaID = aID0(binID);
    
    %Return in the original shape
    subaID_(aID_~=0 & isfinite(aID_)) = subaID;
    subaID = subaID_;
    
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% LINE: RETRIEVE SPECIFIC GENERATIONS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    function [aID,familyMap,newCellFlag,gen0IDs,genMap,A] = getGenerationsT(obj,K,timePt,firstGen,lastGen,useFilt,catEndPts)
    %[aID,familyMap,newCellFlag,gen0IDs,genMap,A] = getGenerationsT(K,timePt,firstGen,lastGen,useFilt).
    %Retrieves the indices of the firstGen through lastGen generations of cells
    %from a given timePt (in minutes) in the Kth line. Note that the zeroth 
    %generation is the generation that is born before timePt and divides after 
    %timePt. useFilt applies the existing data filter to remove cells and
    %any progeny from the lineage. Note that aID{i}{j} is a vector
    %containing the indices in A corresponding to the lineage denoted by
    %familyMap{i}(:,j). All aID{i} cells and familyMap{i} columns have the
    %same zeroth generation cellID. gen0IDs are the cellIDs of the zeroth
    %generation.genMap is the same shape as familyMap, with the generation
    %of each cellID in familyMap indicated.
    
    %Load data from the Kth line
    A = load(fullfile(obj.dirList.flatDir,obj.lineFiles{K}));
    A = A.A;
    
    %Lookup line information if filtering required
    if useFilt
        fov = obj.lineFiles{K};
        line = str2double(fov(4:7));
        fov = str2double(fov(1:2));
    end
    
    %Find the zeroth generation
    %Check the first and last time points of each cell in line K
    diffID = diff(A(:,obj.a('cellID')));
    firstID = [true; logical(diffID)];
    lastID = [logical(diffID); true];
    cellList1 = A(firstID,obj.a('cellID'));
    cellList2 = A(lastID,obj.a('cellID'));
    
    bornBefore = cellList1(A(firstID,obj.a('t')) <= timePt);
    dividesAfter = cellList2(A(lastID,obj.a('t')) >= timePt);
    
    gen0IDs = bornBefore(ismember(bornBefore,dividesAfter));
    nGen0 = numel(gen0IDs);
    
    if useFilt
        %Remove zeroth gen cells if filter requested
        rmCell = false(1,nGen0);        
        bID = obj.getBIndex(fov,line,gen0IDs);
        rmCell(not(bID)) = true; %cells that are not in B (these should not exist though)
        rmCell(logical(bID)) = not(obj.filterIDFlat(bID(bID~=0)));
        
        gen0IDs = gen0IDs(~rmCell);
        nGen0 = numel(gen0IDs);
    end    
    
    aID = cell(1,nGen0);
    familyMap = cell(1,nGen0);
    newCellFlag = cell(1,nGen0);
    genMap = cell(1,nGen0);
    
    AcellID = unique(A(:,[obj.a('cellID') obj.a('parentID')]),'rows');
    for j = 1:nGen0
        %Trace the lineage forward
        fMj = obj.getLineageForward(AcellID,gen0IDs(j),lastGen+1);
        
        %Trim to the first generation desired (if greater than zeroth)
        if firstGen > 0
            fMj = fMj(1+firstGen:end,:);
        end
        
        %Trace back generations if requested (i.e. firstGen<0). Note that
        %every column of fMj has the same heritage line if firstGen<0 since
        %so far, the first row is constant gen0IDs(j).
        if firstGen < 0
            heritage = obj.getLineageReverse(A(:,[obj.a('cellID') obj.a('parentID')]),gen0IDs(j),abs(firstGen));
            if length(heritage)>1
                %Reverse lineage exists
                heritageStart = min(abs(firstGen),length(heritage)-1);
                heritage = heritage(end-heritageStart:end-1);
                
                %Check for heritage cells that need to be filtered out
                if useFilt
                    rmCell = false(numel(heritage),1);
                    bID = obj.getBIndex(fov,line,heritage);
                    rmCell(not(bID)) = true;
                    rmCell(logical(bID)) = not(obj.filterIDFlat(bID(bID~=0)));
                    heritage(rmCell) = 0;
                    heritageStart = find(heritage==0,1,'last');
                    
                    if numel(heritageStart)==0
                        %No heritage cells filtered out
                        heritage = heritage*ones(1,size(fMj,2));
                    elseif heritageStart~=length(heritage) 
                        %At least one generation < 0 passes filter
                        heritage = heritage(heritageStart+1:end);
                        heritage = heritage*ones(1,size(fMj,2));
                    else
                        %-1 gen cell and therefore all previous generations must be filtered out
                        heritage = [];
                    end
                end

                %Concatenate
                fMj = [heritage; fMj];
            end
        else
            heritage = [];
        end
        
        cellIDs = unique(fMj(fMj~=0));
                
        %Filter rest of cells if necessary
        if useFilt
            %Find cells to remove
            rmCell = false(numel(cellIDs),1);            
            bID = obj.getBIndex(fov,line,cellIDs);
            rmCell(not(bID)) = true;
            rmCell(logical(bID)) = not(obj.filterIDFlat(bID(bID~=0)));
            
            %Remove cells
            fMj(ismember(fMj,cellIDs(rmCell))) = 0;
            
            %Remove progeny
            fMj = fMj.*logical(cumprod(fMj));
        end
        
        %Truncate lineage if necessary
        cutoff = find(sum(fMj,2)==0,1,'first');
        if sum(cutoff)~=0
            fMj = fMj(1:cutoff-1,:);
        end
        
        %Remove redundant lines: identical lineages
        fMj = unique(fMj','rows')';
        
        %Remove redundant lines: shortened lineages
        for m = 1:size(fMj,2)
            lastChild = find(fMj(:,m)~=0,1,'last');
            if isempty(lastChild)
                fMj(:,m) = 0;
                continue
            end
            
            [~,c] = find(fMj==fMj(lastChild,m));
            c = c(c~=m);
            if any(c)
                %There is a lineage that tracks at least one generation
                %past the mth lineage
                fMj(:,m) = 0;
            end
        end
        
        fMj = fMj(:,sum(fMj,1)~=0);
        
        gMj = (-size(heritage,1):size(fMj,1)-size(heritage,1)-1)';
                
        %Save whatever of the family map is left
        if ~isempty(fMj)            
            familyMap{j} = fMj;
            genMap{j} = gMj;
            
            %Find the relevant indices in A
            aIDj = containers.Map('keytype','double','valuetype','any');
            for m = 1:numel(cellIDs)
                aIDj(cellIDs(m)) = find(A(:,obj.a('cellID'))==cellIDs(m));
            end
            aIDj(0) = [];

            %Build trajectory of indices in A for each line
            aID{j} = cell(1,size(familyMap{j},2));
            newCellFlag{j} = cell(1,size(familyMap{j},2));
            for m = 1:size(familyMap{j},2)
                aIDm = [];
                for p = 1:size(familyMap{j},1)
                    aIDm = [aIDm; aIDj(familyMap{j}(p,m))]; %#ok<AGROW>
                end
                aID{j}{m} = aIDm;

                newCellFlag{j}{m} = A(aIDm,obj.a('cellID')).*[true; diff(A(aIDm,obj.a('cellID')))~=0];
            end   
        end
    end
    
    %Concatenate endpoints if necessary
    if catEndPts
        [aID,familyMap,newCellFlag] = obj.generationsEndPtCat(A,aID,familyMap,newCellFlag);
    end
    
    end
    
    
    
    
    
    function [aID,familyMap,newCellFlag,gen0IDs,genMap,genaID,A] = getGenerationsT2(obj,K,timePt,firstGen,lastGen,useFilt,catEndPts)
    %[aID,familyMap,newCellFlag,gen0IDs,genMap,genaID,A] = getGenerationsT2(obj,K,timePt,firstGen,lastGen,useFilt,catEndPts)

    %Load A matrix
    A = obj.loadLineK(K);

    %Construct cell array to retrieve indices in A
    [cell2aid,AcellID] = obj.cell2aid(A,K,useFilt);

    %Find zeroth generation
    %Check the first and last time points of each cell in line K
    newCell0 = find([1; diff(AcellID(:,1))]~=0);
    newCell1 = find([diff(AcellID(:,1)); 1]~=0);
    cellList1 = AcellID(newCell0,1);
    cellList2 = AcellID(newCell1,1);

    bornBefore = cellList1(A(newCell0,obj.a('t')) <= timePt);
    dividesAfter = cellList2(A(newCell1,obj.a('t')) >= timePt);

    gen0IDs = bornBefore(ismember(bornBefore,dividesAfter));
    gen0IDs(gen0IDs==0) = [];

    nGen0 = numel(gen0IDs);

    aID = cell(1,nGen0);
    familyMap = cell(1,nGen0);
    newCellFlag = cell(1,nGen0);
    genMap = cell(1,nGen0);
    genaID = cell(1,nGen0);

    if nGen0==0
        return
    end

    for j = 1:nGen0
        %Trace lineage forward
        if lastGen >= 0
            fMj = obj.getLineageForward(AcellID,gen0IDs(j),lastGen+1);
        else
            fMj = gen0IDs(j);
        end
            
        %Trim end if not enough generations were tracked forward
        genEnd = sum(fMj,2);
        genEnd = find(genEnd==0,1);
        if ~isempty(genEnd)
            fMj = fMj(1:genEnd-1,:);
        end

        %Trim to the first generation desired, if greater than zeroth
        if firstGen > 0
            fMj = fMj(1+firstGen:end,:);
        end

        %Trace back generations, if less than zeroth
        if firstGen < 0
            heritage = obj.getLineageReverse(AcellID,gen0IDs(j),abs(firstGen));
            if length(heritage)>1
                %Reverse lineage exists
                heritage = heritage(1:end-1);
                heritage = heritage*ones(1,size(fMj,2));
            else
                %No previous generations tracked (or, if applied, that passed
                %the filter)
                heritage = [];
            end
            fMj = [heritage; fMj]; %#ok<AGROW>
        else
            heritage = [];
        end

        %Store the family map and generational map
        if firstGen <= 0
            gMj = ((-size(heritage,1):size(fMj,1)-size(heritage,1)-1)');
        else
            gMj = ((0:size(fMj,1)-1) + firstGen)';
        end
        
        fMj = fMj(gMj<=lastGen,:);
        gMj = gMj(gMj<=lastGen,:);
        
        familyMap{j} = fMj;    
        genMap{j} = gMj;
        
        aID{j} = obj.familyMap2aID(fMj,cell2aid);

        ncf = aID{j};
        
        if ~isempty(ncf)
            ncf(ncf~=0) = AcellID(ncf(ncf~=0),1);
            ncf = ncf.*[ones(1,size(ncf,2)); diff(ncf,1,1)];
            
            gaid = logical(ncf);
            gaid = cumsum(gaid) + genMap{j}(1) - 1;
            gaid(aID{j}==0) = NaN;            
        else
            ncf = [];
            gaid = [];
        end
        
        newCellFlag{j} = ncf;
        genaID{j} = gaid;    
    end
    
    %Remove empty cells
    keepMe = ~cellfun(@isempty,aID);
    
    gen0IDs = gen0IDs(keepMe);
    aID = aID(keepMe);
    familyMap = familyMap(keepMe);
    newCellFlag = newCellFlag(keepMe);
    genMap = genMap(keepMe);
    genaID = genaID(keepMe);
    
    if catEndPts
        [aID,newCellFlag,genaID] = obj.catEndPts(K,A,aID,newCellFlag,genaID);
    end

    end
    
    
    
    
    
    
    
    
    
    function k = getLineageForward(obj,AcellID,cellID,nGen)
    %k = getLineageForward(AcellID,cellID,nGen)
    %Gets a matrix of cell identifiers for lineages (one column for each 
    %line) proceeding from the cell corresponding to cellID using AcellID, 
    %the subset of A (or B) where the first column contains cellIDs and 
    %the second column the cellID of each cell's parent (if it exists). 
    %Zeros inserted if a lineage does not have a tracked daughter.
    
    AcellID = unique(AcellID,'rows');
    
    lastCellList = []; %List of the nGen-th cells
    
    cellList = cell(1,nGen);
    cellList{1} = cellID;
    
    for j = 2:nGen
        parentsj = cellList{j-1};
        
        childrenj = ismember(AcellID(:,2),parentsj);
        childrenj = AcellID(childrenj,1);
        cellList{j} = childrenj;
        
        lastCellsj = ~ismember(parentsj,AcellID(:,2));
        lastCellList = [lastCellList; parentsj(lastCellsj)]; %#ok<AGROW>
    end
    
    lastCellList = [lastCellList; cellList{end}];
    lastCellList = unique(lastCellList);
    
    %To rebuild lineage, trace backwards from the last cells, using
    %obj.getLineage, with the killFlag set to cellID
    k = zeros(nGen,numel(lastCellList));
    k(1,:) = cellID;
    for j = 1:numel(lastCellList);
        kj = obj.getLineage(AcellID,lastCellList(j),cellID);
        k(2:1+numel(kj),j) = kj(end:-1:1);
    end
    
    k = sortrows(k');
    k = k';
    
    end
    
    
    function [aID,familyMap,newCellFlag] = generationsEndPtCat(obj,A,aID,familyMap,newCellFlag)
    %Concatenates the first and last points of each family line with the
    %last and first points of the corresponding parent and daughters. If
    %there are two daughters, one is chosen at random. Note that no
    %concatenation is performed on the familyMap structure. 
    %N.B. this function is almost identical to dataMP2.buildTreeEndPtCat,
    %altered to handle the shape of generational lineages.
    
    for k = 1:numel(aID)
        for j = 1:numel(aID{k})
            hasParent = A(aID{k}{j}(1),obj.a('parentID'));
            if hasParent~=0
                aIDparent = find(A(:,obj.a('cellID'))==hasParent,1,'last');
                aID{k}{j} = [aIDparent; aID{k}{j}];
                newCellFlag{k}{j} = [hasParent; newCellFlag{k}{j}];
            end
            hasChild = A(:,obj.a('parentID'))==A(aID{k}{j}(end),obj.a('cellID'));
            hasChild = A(hasChild,obj.a('cellID'));
            if numel(hasChild)>0
                rC = randi(numel(hasChild),1,1);
                hasChild = hasChild(rC);
                aIDchild = find(A(:,obj.a('cellID'))==hasChild,1,'first');
                aID{k}{j} = [aID{k}{j}; aIDchild];                
                newCellFlag{k}{j} = [newCellFlag{k}{j}; hasChild];
            end
        end    
    end
    
    end
    
    
    
    function [cell2aid,AcellIDFiltered,AcellIDRaw] = cell2aid(obj,A,K,useFilt)
    %[cell2aid,AcellIDFiltered,AcellIDRaw] = cell2aid(obj,A,K,useFilt)
    %Constructs cell2aid, which is a cell array where the kth cell contains
    %all the indices in A that correspond to cellID = k. If useFilt is set
    %to true (1), cells that do not pass the filter will have empty cells
    %in cell2aid. Note that K is the Kth line in obj.lineFiles.
    %AcellIDFiltered    A(:,[cellID parentID]), with zeros for cells that
    %                   do not pass the current filter
    %
    %AcellIDRaw         A(:,[cellID parentID]), unfiltered
    
    cellIDs0 = unique(A(:,obj.a('cellID')));
    
    if useFilt
        [fov,line] = obj.fovLineK(K);
        bID = obj.getBIndex(fov,line,cellIDs0);
        
        useMe = obj.filterIDFlat(bID);
        skipMe = cellIDs0(~useMe);
        cellIDs = cellIDs0(useMe);
    else
        skipMe = [];
        cellIDs = cellIDs0;
    end
    
    %Remove filtered cells from consideration
    AcellID0 = A(:,[obj.a('cellID') obj.a('parentID')]);
    AcellID = AcellID0;
    skipA = ismember(AcellID(:,1),skipMe);
    AcellID(skipA,:) = 0;
    skipA = ismember(AcellID(:,2),skipMe);
    AcellID(skipA,2) = 0;
    
    %Bin indices by the cellIDs
    aIndex = (1:size(AcellID,1))';
    newCell0 = find([1; diff(AcellID(:,1))]~=0);
    newCell1 = find([diff(AcellID(:,1)); 1]~=0);
    idLength = newCell1 - newCell0 + 1;
    
    aIndex = mat2cell(aIndex,idLength,1);
    
    z = [1; cumsum(idLength)+1];
    z = z(1:end-1);
    binID = AcellID(z,1);
    
    aIndex = aIndex(binID~=0);
    binID = binID(binID~=0);
        
    cell2aid = cell(max(cellIDs),1);
    cell2aid(binID) = aIndex; %cell array containing indices of matrix A, indexed by cellID
    
    AcellIDFiltered = AcellID;
    AcellIDRaw = AcellID0;
    
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% CALCULATING CROSS-CORRELATIONS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    function [meanX,sdX,meanLAGS,nc,X,LAGS,dtMean] = xcorrDataMP2(obj,prop1,prop2)
    
    useFilt = 1;
    catEndPts = 0;
    
    X = cell(1e4,1);
    LAGS = cell(1e4,1);
    
    c = 1;
    
    dtMean = 0;
    
    for k = 1:numel(obj.lineFiles)
        [aID,~,~,A] = obj.buildTree2(k,'family',useFilt,catEndPts);

        for j = 1:numel(aID)
            uniquePts = [ones(size(aID{j},1),1) diff(aID{j},1,2)~=0].*(aID{j}~=0);
            uniquePts = logical(uniquePts);
            
            %Find last common ancestor to ignore parallel trajectories
            [lca,~] = find(diff([uniquePts; zeros(1,size(uniquePts,2))],1,1)==1);
            lca = [sum(uniquePts(:)); lca]; %#ok<AGROW>
                        
            colID = repmat(1:size(aID{j},2),size(aID{j},1),1);
            colID = colID(uniquePts)';
            colID = colID(:);
            
            %Lag vector by index
            frameLag = repmat((1:size(aID{j},1))',1,size(aID{j},2));
            frameLag = frameLag(uniquePts);
            frameLag = frameLag(:);
            frameLag = frameLag*ones(1,numel(frameLag)) - ones(numel(frameLag),1)*frameLag';
                        
            uniquePts = aID{j}(uniquePts);
            uniquePts = uniquePts(:);
            
            %Ignore parallel trajectories for lag consideration
            useUpToMe = sub2ind([numel(uniquePts) numel(uniquePts)],lca(colID),(1:numel(uniquePts))');
            useMe = true(numel(uniquePts));
            useMe(useUpToMe) = false;
            useMe = cumprod(useMe);
            useMe(useUpToMe) = true;
            
            %Build relevant diagonal indices to include self-contained lags
            fsp = find([1; diff(colID)~=0]);
            lsp = find([diff(colID)~=0; 1]);
            
            for d = 1:numel(fsp)
                useMe(fsp(d):lsp(d),fsp(d):lsp(d)) = true;
            end
            
            useMe = logical(useMe);
            
            %Get relevant properties
            P1 = obj.getPropAIndex(uniquePts,A,prop1,k);
            P2 = obj.getPropAIndex(uniquePts,A,prop2,k);
            
            P10 = ones(numel(P1),1)*P1';
            P20 = P2*ones(1,numel(P2));
            
            P10 = P10(useMe);
            P20 = P20(useMe);
            
            %Calculate the mean and standard deviation for different lags
            lags = frameLag(useMe);
                        
            lagBins = sort(unique(lags));
            lagBins = [lagBins(:); lagBins(end)+1];
            [~,~,binID] = histcounts(lags,lagBins);
            
            mu1 = accumarray(binID,P10,[numel(lagBins)-1,1],@nanmean);
            mu2 = accumarray(binID,P20,[numel(lagBins)-1,1],@nanmean);
            
            sd1 = accumarray(binID,P10.^2,[numel(lagBins)-1,1],@nanmean);
            sd2 = accumarray(binID,P20.^2,[numel(lagBins)-1,1],@nanmean);
                                                    
            sd1 = real(sqrt(sd1 - mu1.^2));
            sd2 = real(sqrt(sd2 - mu2.^2));
            
            mu1 = mu1(binID);
            mu2 = mu2(binID);
                                   
            sd1 = sd1(binID);
            sd2 = sd2(binID);
            
            mu10 = zeros(size(frameLag));
            mu10(useMe) = mu1;
            mu1 = mu10;
            
            mu20 = zeros(size(frameLag));
            mu20(useMe) = mu2;
            mu2 = mu20;
            
            sd10 = ones(size(frameLag));
            sd10(useMe) = sd1;
            sd1 = sd10;
            
            sd20 = ones(size(frameLag));
            sd20(useMe) = sd2;
            sd2 = sd20;
                        
            P12 = (ones(numel(P1),1)*P1' - mu1).*(P2*ones(1,numel(P2)) - mu2)./(sd1.*sd2);
            x = P12(useMe);
                        
            %Each line has its own cross-correlation
            lags = lags(~isnan(x) & isfinite(x));
            x = x(~isnan(x) & isfinite(x));
            
            if isempty(lags)
                continue
            end
            
            lagBins = sort(unique(lags));
            lagBins = [lagBins(:); lagBins(end)+1];
            [~,~,binID] = histcounts(lags,lagBins);
            
            z = zeros(numel(lagBins)-1,1);
            z = accumarray(binID,x(:),size(z),@nanmean,0);
                        
            longestCol = logical(aID{j});
            longestCol = sum(longestCol);
            longestCol = find(longestCol==max(longestCol),1);
            
            dt = mean(diff(obj.getPropAIndex(aID{j}(:,longestCol),A,'t')));
            dtMean = dtMean + dt;
                      
            X{c} = z;
            LAGS{c} = lagBins(1:end-1);
            c = c + 1;
        end
    end
    
    dtMean = dtMean/(c-1);
    
    X = X(1:c-1);
    LAGS = LAGS(1:c-1);
    
    x = cell2mat(X);
    lags = cell2mat(LAGS);
        
    meanLAGS = min(lags):max(lags)+1;
    
    [nc,meanLAGS,binID] = histcounts(lags,meanLAGS);
    meanX = zeros(numel(meanLAGS)-1,1);
    meanX = accumarray(binID(binID~=0),x(binID~=0),size(meanX),@nanmean);
    sdX = accumarray(binID(binID~=0),x(binID~=0).^2,size(meanX),@nanmean);
        
    sdX = sdX - meanX.^2;
    sdX = sqrt(sdX);
    
    meanLAGS = meanLAGS(1:end-1);
    
    nc = nc(:);
    
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% DISPLAY RESULTS OF CELL FILTRATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function T = dispCellStats(obj)
    %T = dispCellStats(obj)
    %Displays the cell statistics that result from the current lab state,
    %i.e. with the current directory and data filters in place. Output is a
    %a table containing the cell statistics
    
    %First filters to display
    filterOrder0 = {'nofilt','wholeCC','symmDivTolP'};
        
    %Load data and retrieve cell stats for each replicate
    T = struct;
    
    %Check which filters are being used
    c = obj.filtersInUse;
    flds = fields(c);
    filterOrder = [filterOrder0(:); sort(flds(~ismember(flds,filterOrder0)))];

    %Fill empty filter fields with NaN
    Tflds = fields(T);
    fldDone = false(1,numel(Tflds));

    %Retrieve cell numbers at each filter point
    F = true(size(obj.filterIDFlat));
    T.none = size(obj.flatData.B,1);
    for f = 2:numel(filterOrder)
        if ~any(ismember(filterOrder{f},fields(T)))
            T.(filterOrder{f}) = [];
        end

        if any(ismember(flds,filterOrder{f}))
            F = F & obj.getFlatFilter(filterOrder{f},c.(filterOrder{f}));
            T.(filterOrder{f}) = [T.(filterOrder{f}); sum(F)]; 
        else
            T.(filterOrder{f}) = [T.(filterOrder{f}); NaN];
        end

        if any(ismember(filterOrder{f},Tflds))
            fldDone(f) = true;
        end
    end

    for f = find(~fldDone)
        T.(Tflds{f}) = [T.(Tflds{f}); NaN];
    end
    
    T = struct2table(T,'rownames',{obj.dirList.baseDir});
    
    %Display table of filter stats
    fprintf('\n\t---NCELLS AFTER FILTRATION---\n')
    disp(T)
    
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% LOOKUP INFORMATION FOR HANDLING MULTIPLE DATASETS
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    function T = getShiftTime(obj,media1,media2)
    %T = GETSHIFTTIME(obj,media1,media2).
    %Gets the shift time (in minutes) from media1 to media2.
    
    M = strfind(obj.media,sprintf('%s,%s',media1,media2));
    C = strfind(obj.media,',');
    MC = zeros(1,numel(obj.media));
    MC(C) = 1;
    MC(1) = 1;
    MC = cumsum(MC); 
    T.tInt = zeros(numel(M),1);
    T.tShift = zeros(1,numel(M));
    for m = 1:numel(M)
        T.tShift(m) = obj.switchTime(MC(M(m)));
        if MC(M(m))==1
            T.tInt(m,1) = 0;
        else
            T.tInt(m,1) = obj.switchTime(MC(M(m))-1);
        end
        if MC(M(m))==numel(obj.switchTime)
            T.tInt(m,2) = max(obj.flatData.B(:,obj.b('t0'))+obj.flatData.B(:,obj.b('divT')))+1e-9;
        else
            T.tInt(m,2) = obj.switchTime(MC(M(m))+1);
        end
    end
    end
    
    function T = getMediaTime(obj,media)
    %T = GETMEDIATIME(media)
    %Retrieves the time interval(s) corresponding to a given media.
    
    M = strfind(obj.media,media);
    C = strfind(obj.media,',');
    MC = zeros(1,numel(obj.media));
    MC(C) = 1;
    MC(1) = 1;
    MC = cumsum(MC);
    T = zeros(numel(M),2);
    for m = 1:numel(M)
        if MC(M(m))==1
            T(m,1) = 0;
        else
            T(m,1) = obj.switchTime(MC(M(m))-1);
        end
        if MC(M(m))>numel(obj.switchTime)
            T(m,2) = max(obj.flatData.B(:,obj.b('t0'))+obj.flatData.B(:,obj.b('divT')))+1e-9;
        else
            T(m,2) = obj.switchTime(MC(M(m)));
        end            
    end
    
    
    end
    
    function media = getOvernight(obj)
    %media = GETOVERNIGHT(OBJ)
    yy = obj.dirList.baseDir(1:4);
    yy = str2double(yy);
    if yy < 2016
        media = 'LB';
    else
        media = 'M9S';
    end
    
    end
    
    
    
    
    function [tf,SZ] = shapeConstant(obj,prop,AorB,testSize)
    %[tf,SZ] = shapeConstantA(prop,AorB,testSize)
    %Checks if, after a property is retrieved by getPropAIndex or
    %getPropBIndex (AorB = 'A' or 'B') returns output the same size as the
    %given input indices. The output size is given in SZ. A test index size
    %testSize can be input to see the corresponding output size.
        
    if nargin < 4
        %Dummy indices to test
        idx = [1 0; NaN 1; 1 1];
    else
        idx = ones(testSize);
        if numel(idx)>3
            %Insert some faulty indices
            idx(1) = 0;
            idx(3) = NaN;
        end
    end
            
    
    switch AorB
        case {'A','a'}
            P = obj.getPropAIndex(idx,obj.loadLineK(1),prop,1);
        case {'B','b'}
            P = obj.getPropBIndex(idx,prop);
    end
    
    SZ = size(P);
    tf = all(size(P)==size(idx));
    
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %                   WRITE DATA TO TXT 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [subdir,Ntot] = writedata(obj,propsA,propsB,destdir)
    %[subdir,Ntot] = WRITEDATA(OBJ,PROPSA,PROPSB,DESTDIR)
    %Writes the properties in propsA (cell lines) and propsB (flat data) to 
    %TXT in destdir. If no input is given for destdir, uigetdir is 
    %triggered. Note that by default propsA and propsB will include fov,
    %line, cellID, and parentID
    if nargin<4, destdir = uigetdir; end
    
    groupbycell = false;    %parameter for single cell files
    catendpoints = false;   %------------"------------------
    
    props0 = {'fov','line','cellID','parentID'}; %Default properties to include
    
    %%%%%%%%%%%%%%%%% PARSE PROPSA and PROPSB %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~iscell(propsA), propsA = {propsA}; end
    if ~iscell(propsB), propsB = {propsB}; end
    
    propsA = unique([props0,propsA],'stable');
    propsB = unique([props0,propsB],'stable');
    
    %%%%%%%%%%%%%%%%%%% MAKE DIRECTORIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    subdir0 = fullfile(destdir,obj.strain);
    if ~exist(subdir0,'dir'), mkdir(subdir0); end
    subdir = fullfile(subdir0,obj.dirList.baseDir);
    if ~exist(subdir,'dir'), mkdir(subdir); end
    linedir = fullfile(subdir,'line_data');
    if ~exist(linedir,'dir'), mkdir(linedir); end
    
    %%%%%%%%%%%%%%%%%%%%% WRITE FLAT FILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Pf = obj.flat_cellProps(propsB);
    Pf = struct2cell(Pf);
    Pf = cell2mat(reshape(Pf,1,numel(Pf)));
    Ntot = size(Pf,1);
    
    flat_name = fullfile(subdir,[obj.dirList.baseDir, '_flat.csv']);
    fid_flat = fopen(flat_name,'w');                            %Open file
    fprintf(fid_flat,headerformat(numel(propsB)),propsB{:});    %Write header
    fprintf(fid_flat,formatstr(propsB),Pf');                    %Write data
    fclose(fid_flat);                                           %Close file
    
    %%%%%%%%%%%%%%%%%%%%% WRITE LINE FILES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Pc = obj.singleCC_cellProps(propsA,groupbycell,catendpoints);
    Pc = struct2cell(Pc);
    Pc = cell2mat(reshape(Pc,1,numel(Pc)));
    
    %Separate files into tracked lines
    fov_col = find(strcmp('fov',propsA));
    line_col = find(strcmp('line',propsA));
    fovline = unique(Pc(:,[fov_col,line_col]),'rows');
    for k = 1:size(fovline,1)
        %Find relevant cells for the k-th file
        idk = (Pc(:,fov_col)==fovline(k,1)) & (Pc(:,line_col)==fovline(k,2));
        Pck = Pc(idk,:);
        
        %Write to file
        line_name = fullfile(linedir,sprintf('fov%2.2d_line%4.4d.csv',fovline(k,1),fovline(k,2)));
        fid_line = fopen(line_name,'w');                            %Open file
        fprintf(fid_line,headerformat(numel(propsA)),propsA{:});    %Write header
        fprintf(fid_line,formatstr(propsA),Pck');                   %Write data
        fclose(fid_line);                                           %Close file
        
        %Print output to keep track of writing
        if round(k/50)==k/50
            fprintf('\t%d of %d lines written.\n',k,size(fovline,1));
        end
    end

    %%%%%%%%%%%%%%%%%% FUNCTIONS TO FORMAT STRINGS %%%%%%%%%%%%%%%%%%%%%%%%
    function str = formatstr(props)
        %str = formatstr(props)
        %Outputs the string format for fprintf according to the properties
        %being written
        if ~iscell(props); props = {props}; end
        str = '';
        for pp = 1:numel(props)
            switch props{pp}
                case {'cellID','parentID','fov','line','frame'}
                    str = [str,'%d,']; %#ok<AGROW>
                otherwise
                    str = [str,'%f,']; %#ok<AGROW>
            end
        end
        str = [str(1:end-1),'\n'];
    end
       
    function str = headerformat(n)
    %str = headerformat(n)
    %Output string format for n headers
    str = repmat('%s,',1,n);
    str = [str(1:end-1),'\n'];
    end
    
    %%%%%%%%%%%%%%%%%%%%%%% WRITE README %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    readme_name = fullfile(subdir,'readme.txt');
    fid_readme = fopen(readme_name,'w');
    
    fprintf(fid_readme,'REPLICATE: %s\n',obj.dirList.baseDir);
    fprintf(fid_readme,'STRAIN: %s\n',obj.strain);
    fprintf(fid_readme,'ORIENTATION: %s\n',obj.orientation);
    fprintf(fid_readme,'OVERNIGHT MEDIA: %s\n',obj.getOvernight);
    shiftT = num2str(obj.switchTime,'%d,');
    shiftT = shiftT(1:end-1);
    fprintf(fid_readme,'SWITCH TIME (min): %s\n',shiftT);
    fprintf(fid_readme,'MEDIA: %s\n',obj.media);
    fprintf(fid_readme,'TIME SHIFTED (min): %d\n',obj.currShift);
    fprintf(fid_readme,'COMMENTS: %s\n\n',obj.comments);
    fprintf(fid_readme,'FILTERS USED\n');
    
    fiu = obj.filtersInUse;
    fiu_fields = fields(fiu);
    for ff = 1:numel(fiu_fields)
        if numel(fiu.(fiu_fields{ff}))==1
            fprintf(fid_readme,'\t%s: %d\n',fiu_fields{ff},fiu.(fiu_fields{ff}));
        else
            fprintf(fid_readme,'\t%s: [%f,%f]\n',fiu_fields{ff},fiu.(fiu_fields{ff}));
        end
    end
    
    fprintf(fid_readme,'\nSMOOTHING METHOD USED: %d\n',obj.smoothingMethod.dosmooth);
    sfa = fields(obj.smoothingMethod.A);
    sfb = fields(obj.smoothingMethod.B);
    for ff = 1:numel(sfa)
        fprintf(fid_readme,'\t%s: %s\n',sfa{ff},obj.smoothingMethod.A.(sfa{ff}));
    end
    for ff = 1:numel(sfb)
        fprintf(fid_readme,'\t%s: %s\n',sfb{ff},obj.smoothingMethod.B.(sfb{ff}));
    end
    
    fprintf(fid_readme,'\nN_TOTAL: %d',size(Pc,1));
    
    %Close the file
    fclose(fid_readme);
        
    end
    
    
    
    function obj = expfitProp(obj,AB,prop,nFit)
    %obj = expfitProp(obj,AB,prop,nFit)
    %Applies smoothing filter to retrieved datasets. prop is variable to be
    %smoothed, using nFit points for the fitting window. prop can be
    %a string or a cell array of strings, and nFit must either be a 
    %constant or contain the same size as prop. AB is the flag to indicate
    %a property calculated either on cell lines (A) or flat data (B), and
    %must either be a single character or the same size as prop.
    
    if ~iscell(prop), prop = {prop}; end
    
    if numel(AB)==1; AB = {AB}; AB = repmat(AB,1,numel(prop)); end
    if numel(nFit)==1; nFit = repmat(nFit,1,numel(prop)); end
    
    %Turn on smoothing filter
    obj.smoothingMethod.dosmooth = true;
    
    for k = 1:numel(prop)
        obj.smoothingMethod.(AB{k}).(prop{k}) = ...
            sprintf('expfit%s%%s%d_0',AB{k},nFit(k));
    end
    
    end
    
    
    
    end %METHODS: Visible, non-static

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    methods (Static = true)
    
    function k = getLineage(AcellID,cellID,killFlag)
    %k = obj.getLineage(AcellID,cellID,killFlag)
    %Gets the cell identifiers of each parent, grandparent, great-grandparent,
    %etc. for the cell corresponding to cellID using AcellID, the
    %subset of A (or B) where the first column contains cellIDs and the
    %second column the cellID of each cell's parent (if it exists).
    %killFlag is the list of cellIDs that in the case of independent
    %lines, have already been saved, therefore the line should
    %terminate before this point. Set killFlag to Inf otherwise.
        
    k = zeros(100,1);
    k(1) = cellID;
    parentID = AcellID(AcellID(:,1)==cellID,2);
    if isempty(parentID)
        k = k(1);
        return
    elseif killFlag==cellID
        k = [];
        return
    end
    parentID = parentID(1);
    kc = 2;
    while any(parentID) && ~any(parentID==killFlag)
        k(kc) = parentID;
        parentID = AcellID(AcellID(:,1)==parentID,2);
        parentID = parentID(1);
        kc = kc + 1;
    end
    k = k(1:kc-1);

    
    end
    
    
    function [k,shortFlag] = getLineageReverse(AcellID,cellID,nGen)
    %[k,shortFlag] = obj.getLineage(AcellID,cellID,nGen)
    %Gets the cell identifiers of each parent, grandparent, great-grandparent,
    %etc. up to the nGen-th generation for the cell corresponding to 
    %cellID using AcellID, the subset of A (or B) where the first column 
    %contains cellIDs and the second column the cellID of each cell's 
    %parent (if it exists). If lineage tracking falls short of the nGen-th
    %generation, shortFlag will be 1; k will only list tracked cellIDs.
        
    k = zeros(nGen+1,1);
    k(1) = cellID;
    parentID = AcellID(AcellID(:,1)==cellID,2);
    if isempty(parentID)
        k = k(1);
        return
    end
    parentID = parentID(1);
    kc = 2;
    while any(parentID) && kc <= (nGen+1)
        k(kc) = parentID;
        parentID = AcellID(AcellID(:,1)==parentID,2);
        parentID = parentID(1);
        kc = kc + 1;
    end
    k = k(1:kc-1);
    k = k(end:-1:1);
    
    %Check lineage length
    if length(k) < (nGen+1)
        shortFlag = true;
    else
        shortFlag = false;
    end
    
    end
    
    
    
    function aid = familyMap2aID(familyMap,cell2aid)
    %[aid,genaID] = familyMap2aID(familyMap,cell2aid).
    %Uses output from dataMP.cell2aid to reconstruct the indices from an A
    %matrix to correspond to the list of cells in the lineage given by
    %familyMap (a vector).
    
    fMj = familyMap;
    
    %Build index trajectories
    %If only one line, do not need to pad with zeros
    if size(fMj,2)==1
        aid = cell2aid(familyMap);
        aid = cell2mat(aid(:));
        return
    end
    
    %If more than one line, need to reconstruct with padding to account for
    %lineages of different lengths
    aid = cell(size(fMj));
    aid(fMj~=0) = cell2aid(fMj(fMj~=0));
    
    aid = mat2cell(aid,size(aid,1),ones(1,size(aid,2)));
    aid = num2cell(aid);
    aid = cellfun(@(x) cell2mat(x{1}),aid,'uniformout',0);
    
    %Pad tracks of different lengths with zeros
    naid = cellfun(@numel,aid);
    maxnaid = max(naid);
    if ~all(naid==maxnaid)
        aid = [aid; cell(1,numel(aid))];
        padMe = zeros(maxnaid*size(aid,2)-sum(naid),1);
        padMe = mat2cell(padMe,maxnaid-naid,1);
        aid(2,:) = padMe;
        
        aid = mat2cell(aid,size(aid,1),ones(1,size(aid,2)));
        aid = num2cell(aid);
        aid = cellfun(@(x) cell2mat(x{1}),aid,'uniformout',0);
    end
    
    aid = cell2mat(aid);
    
    end
    
    
    function pID = uniqueGen(FMorAID,GMorGAID,genOfInterest)
    %pID = uniqueGen(FMorAID,GMorGAID,genOfInterest)
    %Outputs the index for the property matrix output by e.g.
    %line_cellPropsGenerationT for the generation of interest. Only outputs
    %unique cells.
    
    %Identify cells belonging to the generation of interest
    if size(FMorAID,2)>1 && size(GMorGAID,2)>1
        goi = GMorGAID==genOfInterest;
    else
        goi = repmat(GMorGAID==genOfInterest,1,size(FMorAID,2));
    end
    
    isUnique = [ones(size(FMorAID,1),1) diff(FMorAID,1,2)].*(FMorAID~=0);
    pID = goi & logical(isUnique);
    
    end
    
    
    
    function [otherGenUsed,D] = formatOtherGenUsed(D,useParent,extensiveFlag)
    %[otherGenUsed,D] = formatOtherGenUsed(D,useParent,extensiveFlag)
    %Outputs the matrix otherGenUsed to multiply property matrices obtained
    %during getPropAIndex procedures for e.g. cell cycle smoothing or cross-
    %generational calculations. D is the structure output from getdiffaID.
    %If useParent is set to false, otherGenUsed will be set to false for
    %all points outside of the cell cycle of interest, which overrides the
    %extensiveFlag. If the extensiveFlag is true, daughter points are set
    %to 2 and parents to 1/2.
    
    daID = D.daID;
    parentUsed = D.parentUsed;
    daughterUsed = D.daughterUsed;
    ptRange = D.ptRange;
    
    %If an extensive property: double daughter values and half parent
    if ~useParent
        parentUsed = parentUsed.*1;
        parentUsed(parentUsed==1) = NaN;
        parentUsed(parentUsed==0) = 1;
        parentUsed = cumprod(parentUsed(:,end:-1:1),2);
        parentUsed = parentUsed(:,end:-1:1);
        
        daughterUsed = daughterUsed.*1;
        daughterUsed(daughterUsed==1) = NaN;
        daughterUsed(daughterUsed==0) = 1;
        daughterUsed = cumprod(daughterUsed,2);
        
        otherGenUsed = [parentUsed daughterUsed];
        
        %Recalculate valid point range
        ptRange0 = ptRange;
        [r,c] = find(~isnan(otherGenUsed));
        s = sortrows([r c]);
        ptRange = zeros(size(daID,1),2);
        ptRange(s(end:-1:1,1),1) = s(end:-1:1,2);
        ptRange(s(:,1),2) = s(:,2);
        ptRange(:,1) = max(ptRange(:,1),ptRange0(:,1));
        ptRange(:,2) = min(ptRange(:,2),ptRange0(:,2));
        
        D.ptRange = ptRange;
        
    elseif extensiveFlag
        parentUsed = parentUsed.*0.5;
        parentUsed(parentUsed==0) = 1;
        parentUsed = cumprod(parentUsed(:,end:-1:1),2);
        parentUsed = parentUsed(:,end:-1:1);
        
        daughterUsed = daughterUsed.*2;
        daughterUsed(daughterUsed==0) = 1;
        daughterUsed = cumprod(daughterUsed,2);
        
        otherGenUsed = [parentUsed daughterUsed];
    else
        otherGenUsed = ones(size(daID));
    end
    
    end   
    
    
    function V = volume(L,w,area,mthd)
    %V = volume(L,w,area,mthd)
    %Calculates the volume of a cell using a specific method. Set mthd to
    %either 'hemisph' (cyl + hemispherical caps), 'hemiellip' (cyl + 
    %hemi-ellipsoidal caps), 'cyl' (just a cylinder), or 'area' (returns
    %the same area that was input). area input not required field for 
    %either hemispherical or cylinder calculation, set to []. 

    if nargin < 4
        mthd = 'hemisph';
    end
        
    switch mthd
        case 'hemisph'
            V = 4/3*pi*(w/2).^3 + pi*(w/2).^2.*(L-w);
            
        case 'hemiellip'
            %Ellipsoid parameters calculated assuming that the area is a
            %rectangle and ellipse.
            a = w/2;
            b = (L.*w - area)./(w.*(2 - pi/2));
            V = 4/3*pi*a.^2.*b + pi*a.^2.*(L - 2*b);
            
        case 'cyl'
            V = pi*(w/2).^2.*L;
            
        case 'area'
            V = area;
        
        otherwise
            %Default to hemispherical caps + cylinder assumption
            V = 4/3*pi*(w/2).^3 + pi*(w/2).^2.*(L-w);
    end
    
    end
    
    
    
    function SA = surfacearea(L,w,area,mthd)
    %SA = surfacearea(L,w,area,mthd)
    %Calculates the surface area of a cell using a specific method. Set mthd to
    %either 'hemisph' (cyl + hemispherical caps), 'hemiellip' (cyl + 
    %hemi-ellipsoidal caps), 'cyl' (just a cylinder), or 'area' (returns
    %the same area that was input). area input not required field for 
    %either hemispherical or cylinder calculation, set to []. 
    
    if nargin < 4
        mthd = 'hemisph';
    end
    
    switch mthd
        case 'hemisph'
            SA = 4*pi*(w/2).^2 + pi*w.*(L-w);
                        
        case 'hemiellip'
            %Ellipsoid parameters calculated assuming that the area is a
            %rectangle and ellipse.
            a = w/2;
            b = (L.*w - area)./(w.*(2 - pi/2));
            
            %Different formula if ellipsoid is oblate (a>b) or prolate (a<b)
            SA = zeros(numel(a),1);
            
            ecc = sqrt(1-b.^2./a.^2);
            SA1 = 2*pi*a.^2.*(1+(1-ecc.^2)./ecc.*atanh(ecc));
            
            ecc = sqrt(1-a.^2./b.^2);
            SA2 = 2*pi*a.^2.*(1+b./(a.*ecc).*asin(ecc));
            
            SA(a>b) = SA1(a>b);
            SA(a<b) = SA2(a<b);
            
        case 'cyl'
            SA = pi*w.*L;
                    
        otherwise
            %Default to hemispherical caps + cylinder assumption
            SA = 4*pi*(w/2).^2 + pi*w.*(L-w);
    end
    
    end
    
    
    function w = width(area,L)
    %w = width(area,L)
        w = 2*(L-sqrt(L.^2-(4-pi)*area))./(4-pi);
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    function tf = isExtensive(prop)
    %tf = isExtensive(prop)
    %Checks if a given property is extensive with respect to division, i.e.
    %if it should be halved or doubled when considering a cell trajectory
    %that spans multiple generations.
    
    extensiveProps = {'Itot','area','L','SA'};
    
    if any(strfind(prop,'growthRate'))
        tf = false;
    elseif strcmp(prop(1),'V')        
        %Any calculation of volume is extensive
        tf = true;
    elseif regexp(prop,'^SA')
        %Also for surface area
        tf = true;
    elseif regexp(prop,'^Itot')
        tf = true;        
    elseif any(strcmp(prop,extensiveProps))
        %Matches extensive properties listed above
        tf = true;
    else
        tf = false;
    end
    
    end
    
    
    function [baseProp,mthd,raw] = parseMethod(str)
    %[baseProp,mthd,raw] = parseMethod(str)
    %Separates the string into the base property and the method used to
    %calculate it
    mthdList = {'hemisph','hemiellip','cyl'};
    R = '';
    for r = 1:numel(mthdList)
        R = [R sprintf('%s|',mthdList{r})]; %#ok<AGROW>
    end
    
    baseProp = regexp(str,['(raw|)(\w*?)(' R '$)'],'tokens');
    mthd = baseProp{1}{3};
    raw = baseProp{1}{1};
    baseProp = baseProp{1}{2};
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    function S = getSmoothingParams(prop)
    
    switch prop
        case 'area'
            S = [0.7 0.7];
        case 'L'
            S = [0.78 0.85];
        case {'Vhemisph','V'}
            S = [0.65 0.61];
        case {'SA','SAhemisph'}
            S = [0.72 0.65];    
        otherwise
            S = [];    
    end
    
    end
    
    end
end