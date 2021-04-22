classdef dataLabMP2 < handle
    properties
        dirList
        mediaList
        strainList
        promoterList
        lociList
        dirFilter
        activeDirList
        dataFilter
        smoothingMethod
        strainMap
        RGB
        recordPath
        recordID
        figDir
        figHandles
        figNames
        figUnique
        filterMATpath
        dataDir %directory to save any MATs 
    end
    
    properties (SetAccess = immutable)
        dirFilterFixed
    end
    
    properties (SetAccess = private, Hidden = true)
        recordPathDefault
    end
    
    %%
    methods
        
    function obj = dataLabMP2
        %INITIALIZE LAB OBJECT
        D = allDataBookDirsMP;
        obj.dirList = D.dirList;
        obj.mediaList = D.mediaList;
        obj.strainList = D.strainList;
        obj.promoterList = D.promoterList;
        obj.lociList = D.lociList;
        
        %Filter replicates by strain, loci, overnight, etc.
        obj.dirFilter = struct('strain','','loci','','overnight','','shift','','notReplicate','','notStrain','');
        
        %Fixed filters that are unlikely to be modified
        obj.dirFilterFixed = struct;
        obj.dirFilterFixed.processed  = true;
        obj.dirFilterFixed.notReplicate = {'20160128';'20160129';'20160228'};
        obj.dirFilterFixed.notStrain = {''};
        
        %Initialize filter to select active directories for analysis
        obj.initFilter;
        
        %Smoothing method to use during calculations of e.g. L, area
        obj.initSmoothingMethod;
        
        %Strain map: numeric coding for strains, loci, promoter, etc
        obj.strainMap = obj.initStrainMap;
        obj.RGB = obj.initRGB;
        
        %Directory to save figures, use current directory by default
        obj.initFigDir;
        
        %Flag that ensures all figures have unique names and that
        obj.figUnique = false;
        
        %Path to a record of save instances of the dataLab
        obj = obj.updatePaths;
        obj.recordID = dec2hex(round(now*24*60*60));
        
        %Initialize the record file if it does not already exist
        obj.initRecord;
               
        
    end
    
    function obj = updatePaths(obj)
    %Updates the paths that the lab record and data filters are saved (in
    %case of changed directories or machines. Also updates figure and data
    %save directories if it is currently linked to dropbox on another 
    %machine.
    
    dataHandlingDir = which('dataLabMP2');
    [dataHandlingDir,~,~] = fileparts(dataHandlingDir);
    
    obj.recordPathDefault = fullfile(dataHandlingDir,'labRecord.mat');
    obj.recordPath = obj.recordPathDefault;
    
    %MAT containing filters that are relatively costly to
    %calculate with each initialization of dataMP2
    obj.filterMATpath = fullfile(dataHandlingDir,'dataFilters.mat');
    
    %Update dropbox folders to link correctly
    if any(strfind(obj.figDir,'Dropbox'));
        %Figure and data save directories
        Z = {'figDir', 'dataDir'};
        for k = 1:2
            z = obj.(Z{k});
            if ispc
                z(strfind(z,'/')) = '\';
                z = ['C:\Users\Mia\' z(strfind(z,'Dropbox'):end)];
            else
                z(strfind(z,'\')) = '/';
                z = ['/Users/miapanlilio/' z(strfind(z,'Dropbox'):end)];
            end
            obj.(Z{k}) = z;
        end
    end
    
    end
    
    function obj = initSmoothingMethod(obj)
    %Initialize structure with smoothing method to use for different
    %metrics
    
    obj.smoothingMethod = struct('dosmooth',false,'A',struct,'B',struct);
    
    end
    
    function obj = initFilter(obj)
    %Initialize folder filter
    D = obj.fixedFilter;
    obj.activeDirList = obj.dirList(D);
    
    obj.dataFilter = struct;
    obj.dataFilter.wholeCC = true;
    obj.dataFilter.symmDivTolP = [0.4 0.6];
    
    flds = fields(obj.dirFilter);
    for f = 1:numel(flds)
        obj.dirFilter.(flds{f}) = '';
    end
    
    end
    
    function obj = setDirFilter(obj,filterName,filterVal)
    %OBJ = SETDIRFILTER(OBJ,FILTERNAME,FILTERVAL)
    if ~iscell(filterName), filterName = {filterName}; filterVal = {filterVal}; end
    if ~iscell(filterVal), filterVal = {filterVal}; end
        
    for f = 1:numel(filterName)
        obj.dirFilter.(filterName{f}) = filterVal{f};
    end
    
    D0 = obj.fixedFilter;
    [~,D] = obj.filterDir;
    
    obj.activeDirList = obj.dirList(D0 & D);
    
    end
    
    
    function D = fixedFilter(obj)
    filterNames = fields(obj.dirFilterFixed);
    D = true(numel(obj.dirList),1);
    for f = 1:numel(filterNames)
        [~,D0] = obj.filterDir(filterNames{f},obj.dirFilterFixed.(filterNames{f}));
        D = D & D0;
    end
    
    end
    
    
    function [D,logicalD] = filterDir(obj,filterName,filterVal)
        %Retrieve directories corresponding to a particular strain,
        %overnight media, loci, and/or media condition.
        if nargin==1
            flds = fields(obj.dirFilter);
            useMe = false(numel(flds),1);
            filterVal = cell(numel(flds),1);
            for f = 1:numel(flds)
                if ~isempty(obj.dirFilter.(flds{f}))
                    useMe(f) = true;
                end
                filterVal{f} = obj.dirFilter.(flds{f});
            end
            filterName = flds(useMe);
            filterVal = filterVal(useMe);
        else
            if ~iscell(filterName), filterName = {filterName}; end
            if ~iscell(filterVal), filterVal = {filterVal}; end
        end
        
        D0 = true(numel(obj.dirList),1);
        for j = 1:numel(obj.dirList)
            for k = 1:numel(filterName)
                info = lookupDataMP(obj.dirList{j});
                switch filterName{k}
                    case 'strain'
                        D0(j) = D0(j) && strcmp(info.STRAIN,filterVal{k});
                    case 'loci'
                        D0(j) = D0(j) && strcmp(obj.getLoci(info.STRAIN),filterVal{k});
                    case 'overnight'
                        D0(j) = D0(j) && strcmp(obj.getOvernight(obj.dirList{j}),filterVal{k});
                    case 'shift'
                        D0(j) = D0(j) && any(strfind(info.MEDIA,filterVal{k}));
                    case 'notReplicate'
                        D0(j) = D0(j) && ~any(strcmp(info.FOLDER,[filterVal{k}(:); obj.dirFilterFixed.notReplicate]));
                    case 'notStrain'
                        D0(j) = D0(j) && ~any(strcmp(info.STRAIN,[filterVal{k}; obj.dirFilterFixed.notStrain]));
                    case 'processed'
                        D0(j) = D0(j) && (info.PROCESSED==obj.dirFilterFixed.processed);
                end
            end
        end
        logicalD = D0;
        D = obj.dirList(D0);
    end
    
    
    function obj = storeFigHandle(obj,hf,name)
    %Store figure handles to save later
    if ~iscell(name), name = {name}; end
    
    for hh = 1:numel(hf)
        c = numel(obj.figHandles) + 1;
        obj.figHandles(c) = hf(hh);
        obj.figNames{c} = name{hh};
    end
    
    end
    
    function printFigs(obj)
        
    if obj.figUnique
        rID = dec2hex(round(now*24*60*60));
        
        try
            %Check if there is an exisiting record and check if anything,
            %e.g. filters have changed
            load(obj.recordPath);
            if isKey(labRecord,obj.recordID) %#ok<NODEF>
                %Old record exists
                oldRecord = labRecord(obj.recordID);
                if ~isequal(obj,oldRecord)
                    %If the records do not match, make a new one
                    obj.recordID = rID;
                    labRecord(obj.recordID) = obj; %#ok<NASGU>
                    save(obj.recordPath,'labRecord')
                end
            else
                %Old record does not exist
                rID = obj.recordID;
                labRecord(obj.recordID) = obj; %#ok<NASGU>
                save(obj.recordPath,'labRecord')
            end
            
            str = ['_r' rID];
            
        catch ME
            ME %#ok<NOPRT>
            warning('Lab record ID %s will not be saved.',obj.recordID)
            rID = dec2hex(round(now*24*60*60));
            str = ['_r' rID];
        end
            
    else
        str = '';
    end
    
    for k = 1:numel(obj.figHandles)
        try
            print2pdfMP(obj.figHandles(k),fullfile(obj.figDir,[obj.figNames{k} str]));
            saveas(obj.figHandles(k),fullfile(obj.figDir,[obj.figNames{k} str]),'fig');
        catch ME
            ME %#ok<NOPRT>
        end        
    end
    end
    
    function obj = initFigDir(obj)
        obj.figDir = cd;
        obj.figHandles = [];
        obj.figNames = []; 
        obj.dataDir = cd;
    end
    
    function obj = clearFigHandles(obj)
        obj.figHandles = [];
        obj.figNames = [];
    end
    
    function obj = initRecord(obj)
    %Initialize the lab record file
    if ~exist(obj.recordPath,'file')
        if regexp(obj.recordPath,'.mat$')
            %Valid MAT to write        
            labRecord = containers.Map;
            labRecord(obj.recordID) = obj; %#ok<NASGU>
            save(obj.recordPath,'labRecord');
        else
            %Use the default
            if exist(obj.recordPathDefault,'file')
                %Use the default if it exists
                warning('%s is not a valid lab record. Default %s used.',obj.recordPath,obj.recordPathDefault)
                obj.recordPath = obj.recordPathDefault;
            elseif regexp(obj.recordPathDefault,'.mat$')
                %Create the default if it doesn't exist
                labRecord = containers.Map;
                labRecord(obj.recordID) = obj; %#ok<NASGU>
                save(obj.recordPath,'labRecord');
            else
                %Don't create the lab object
                warning('No valid lab record provided by default path.')
            end 
        end
    end

    end
    
    
    function saveRecord(obj)
        
    load(obj.recordPath);
    if ~isKey(labRecord,obj.recordID) %#ok<NODEF>
        labRecord(obj.recordID) = obj; %#ok<NASGU>
        save(obj.recordPath,'labRecord')
    else
        q = input(sprintf('%s already stored in lab record. Overwrite?',obj.recordID));
        if q
            labRecord(obj.recordID) = obj; %#ok<NASGU>
            save(obj.recordPath,'labRecord')
            fprintf('%s overwritten.\n',obj.recordID);
        else
            fprintf('No new record written.\n');
        end
        clear labRecord            
    end
    
    end
    
    function recordDate(obj)
    %Displays the date that this lab was produced
    disp(datestr(hex2dec(obj.recordID)./(24*60*60)))
    
    end
    
    function loci = getLoci(obj,strain)
        %loci = GETLOCI(obj,infoFile)
        loci = '';
        for k = 1:numel(obj.lociList)
            f = strfind(strain,obj.lociList{k});
            if any(f)
                loci = obj.lociList{k};
                return
            end
        end
    end
    
    function T = dispCellStats(obj)
    %T = dispCellStats(obj)
    %Displays the cell statistics that result from the current lab state,
    %i.e. with the current directory and data filters in place. Output is a
    %a table containing the cell statistics
    
    %Display active filters
    fprintf('\n\t---DIRECTORY FILTERS---\n')
    disp(obj.dirFilter)
    fprintf('\n\t---DATA FILTERS---\n')
    disp(obj.dataFilter)
    
    %First filters to display
    filterOrder0 = {'wholeCC';'symmDivTolP'};
        
    %Load data and retrieve cell stats for each replicate
    T = struct;
    
    for k = 1:numel(obj.activeDirList)
        %Load data
        a = obj.loadDataMP2(obj.activeDirList{k});
        
        %Check which filters are being used
        c = a.filtersInUse;
        flds = fields(c);
        filterOrder = [filterOrder0; flds(~ismember(flds,filterOrder0))];
        
        %Fill empty filter fields with NaN
        Tflds = fields(T);
        fldDone = false(1,numel(Tflds));
        
        %Retrieve cell numbers at each filter point
        F = true(size(a.filterIDFlat));
        for f = 1:numel(filterOrder)
            if ~any(ismember(filterOrder{f},fields(T)))
                T.(filterOrder{f}) = [];
            end
            
            if any(ismember(flds,filterOrder{f}))
                F = F & a.getFlatFilter(filterOrder{f},c.(filterOrder{f}));
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
    end
    
    Tflds = fields(T);
    for f = 1:numel(Tflds)
        T.(Tflds{f}) = [T.(Tflds{f}); nansum(T.(Tflds{f}))];
    end
    
    T = struct2table(T,'rownames',[obj.activeDirList; 'TOTAL']);
    
    %Display table of filter stats
    fprintf('\n\t---NCELLS AFTER FILTRATION---\n')
    disp(T)
    
    end
    
    function T = dispCellStatsWindow(obj,tInt)
    %T = dispCellStatsWindow(obj,tInt)
    %Displays the cell statistics as a result of the current filters, and
    %for a given time interval (tInt), here assuming time is given as 
    %before or after a shift. Includes any cells that are born or divide
    %within the given time window.
    
    T = obj.dispCellStats;
    
    N = zeros(numel(obj.activeDirList)+1,1);
    for k = 1:numel(obj.activeDirList)
        %Load data
        a = obj.loadDataMP2(obj.activeDirList{k});
        
        %Retrieve the final filter vector
        F0 = a.filterIDFlat;
        
        %Check which cells are located in the time interval of interest
        F = a.getFlatFilter('t0Trans',[-Inf tInt(2)]) & a.getFlatFilter('t1Trans',[tInt(1) Inf]) & F0;
        N(k) = sum(F);
    end
    
    N(end) = sum(N);
    
    %Concatenate tables
    N = table(N,'variablenames',{'window'});
    T = [T N];
    
    %Display table of filter stats
    fprintf('\n\t---NCELLS AFTER FILTRATION---\n')
    disp(T);
    
    
    end
    
    function dataMP2obj = loadDataMP2(obj,baseDir)
    %dataMP2obj = loadDataMP2(obj,baseDir)
    %Loads the dataMP2 object. Applies filters and smoothing parameters
    
    %Load data    
    dataMP2obj = dataMP2(getDirsMP_SSD(baseDir));
    
    %Write any smoothing parameters to use
    %NB this field was created in 2017-11: included conditional statements
    %to create these fields for e.g. old record of dataLabMP2
    if isempty(obj.smoothingMethod)
        obj.initSmoothingMethod;
    end
    
    dataMP2obj.smoothingMethod.dosmooth = obj.smoothingMethod.dosmooth;
    
    smthA = fields(obj.smoothingMethod.A);
    for k = 1:numel(smthA)
        dataMP2obj.smoothingMethod.A.(smthA{k}) = obj.smoothingMethod.A.(smthA{k});
    end
    
    smthB = fields(obj.smoothingMethod.B);
    for k = 1:numel(smthB)
        dataMP2obj.B.(smthB{k}) = obj.smoothingMethod.B.(smthB{k});
    end
    
    %Apply any existing filters
    filters = fields(obj.dataFilter);
    
    %Order filters to ensure that wholeCCRobust and omitShift filters are
    %executed last
    if any(strcmp(filters,'wholeCCRobust'))
        filters(strcmp(filters,'wholeCCRobust')) = [];
        filters = [filters; {'wholeCCRobust'}];
    end
    if any(strcmp(filters,'omitShift'))
        filters(strcmp(filters,'omitShift')) = [];
        filters = [filters; {'omitShift'}];
    end
    
    for f = 1:numel(filters)
        if strcmp(filters{f},'steady')
            %Steady state filter
            if isempty(obj.dataFilter.steady)
                continue
            else
                media = obj.dataFilter.steady;
                T = getMediaTime(dataMP2obj,media);
                
                %Grab the first interval if there are multiple
                T = T(1,:);
                
                %Determine whether a shift happened directly before the
                %start of this media interval and if some transition time
                %should be allowed
                transitionTime = 5*60;
                T_ = dataMP2obj.switchTime;
                if any(T_<=T(1,1))
                    T(1,1) = T(1,1) + transitionTime;
                end
                
                dataMP2obj.setFlatFilter({'t0','t1'},{[T(1,1) Inf],[-Inf T(1,2)]});
%                 dataMP2obj.setFlatFilter({'t0','t1'},{[-Inf T(1,2)],[max(0,T(1,1)) T(1,2)]});

                %Shift and steady data are mutually exclusive
                obj.dataFilter.shift = '';
            end
        elseif strcmp(filters{f},'shift')
            %Look at transient from the start of one steady state
            %region to the next
            if isempty(obj.dataFilter.shift)
                continue
            else
                media = regexp(obj.dataFilter.shift,'(\w+),(\w+)','tokens'); 
                Tshift = dataMP2obj.getShiftTime(media{1}{1},media{1}{2});

                if numel(Tshift.tShift)==0
                    %No shift found: all cells are filtered out
                    dataMP2obj.setFlatFilter('t0',[1 -1]) %Impossible filter to pass
                    return
                else
                    %Check the first shift, just in case there were
                    %multiple steps
                    T0 = Tshift.tShift(1);
                    tInt0 = Tshift.tInt(1,:);

                    T1ok = [];%Tsteady1(Tsteady1(:,1) < T0,:);
                    T2ok = [];%Tsteady2(Tsteady2(:,1) > T0,:);

                    tRange = [0 0];
                    if numel(T1ok)~=0
                        tRange(1) = max(tInt0(1),T1ok(1));
                    else
                        %Just take nH hour interval before the shift
                        nH1 = 6;
                        tRange(1) = max(tInt0(1),T0-nH1*60);
                        fprintf('No appropriate steady state found in first media. %d h interval chosen.\n',nH1);
                    end
                    if numel(T2ok)~=0
                        tRange(2) = min(tInt0(2),T2ok(2));
                    else
                        nH2 = 8;
                        tRange(2) = min(tInt0(2),T0+nH2*60);
                        fprintf('No appropriate steady state found in second media. %d h interval chosen.\n',nH2);
                    end

                    %Apply time range filter
                    dataMP2obj.setFlatFilter({'t0','t1'},{[tRange(1) Inf],[-Inf tRange(2)]});
%                     dataMP2obj.setFlatFilter({'t0','t1'},{[-Inf tRange(2)],[tRange(1) Inf]});

                    %Note the shift of interest
                    dataMP2obj.currShift = T0;

                    %Shift and steady data are mutually exclusive
                    obj.dataFilter.steady = '';
                end                     
            end
        elseif strcmp(filters{f},'shiftManual')
            %Look at the shift based on a manually selected region
            %(stored in manualSteadyStates.m)
            if isempty(obj.dataFilter.shiftManual)
                continue
            else
                C = manualSteadyStates(obj.dataFilter.shiftManual);
                tIntManual = C(dataMP2obj.dirList.baseDir);
                dataMP2obj.setFlatFilter({'t0','t1'},{[tIntManual(1) tIntManual(end)],[tIntManual(1) tIntManual(end)]});

                media = regexp(obj.dataFilter.shiftManual,'(\w+),(\w+)','tokens'); 
                Tshift = dataMP2obj.getShiftTime(media{1}{1},media{1}{2});

                dataMP2obj.currShift = Tshift.tShift(1); %Take the first shift
                obj.dataFilter.steady = '';
            end
        else
            dataMP2obj.setFlatFilter(filters{f},obj.dataFilter.(filters{f}));
        end
    end
    
    
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    function hf = replicateOverviewAll(obj)
    hf = zeros(numel(obj.activeDirList),1);
    tic
    for k = 1:numel(obj.activeDirList)
        a = dataMP2(getDirsMP_SSD(obj.activeDirList{k}));
        hf(k) = replicateOverview(a);
        
        obj.storeFigHandle(hf(k),sprintf('overview_%s',obj.activeDirList{k}));
        fprintf('%d of %d done. Time elapsed = %5.5f sec \n',k,numel(hf),toc)
    end
    end
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    function [M,SD,tBins,tRANGE,RHO,hf] = replicateCORRall(obj,prop1,prop2)
    
    tic
    RHO = cell(numel(obj.activeDirList),1);
    tRANGE = cell(numel(obj.activeDirList),1);
    
    ADL = obj.activeDirList;
    parfor k = 1:numel(ADL)
        a = dataMP2(getDirsMP_SSD(ADL{k}));
        a = a.setFlatFilter({'wholeCC','symmDivTolP'},{true,[0.4 0.6]});
        
        [tRange,rho] = replicateCORR(a,prop1,prop2);
        
        RHO{k} = rho;
        tRANGE{k} = tRange(:);
        fprintf('%d of %d done. Time elapsed = %5.5f sec \n',k,numel(obj.activeDirList),toc)
    end
    
    %Get the mean and SD
    RHO = cell2mat(RHO);
    tRANGE = cell2mat(tRANGE);
    tBins = floor(min(tRANGE)/5)*5:5:ceil(max(tRANGE)/5)*5;
    
    [~,~,binID] = histcounts(tRANGE,tBins);
    
    M = accumarray(binID,RHO,[numel(tBins)-1 1],@mean);
    SD = accumarray(binID,RHO.^2,[numel(tBins)-1 1],@mean);
    SD = SD - M.^2;
    SD = sqrt(SD);
    
    hf = figure;
    plot(tBins(1:end-1)+diff(tBins(1:2))/2,M,'k'), hold on
    plot(tBins(1:end-1)+diff(tBins(1:2))/2,M+SD,'color',[1 1 1]*0.7)
    plot(tBins(1:end-1)+diff(tBins(1:2))/2,M-SD,'color',[1 1 1]*0.7)
    hold off
    
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function P = flat_cellProps(obj,props)
    if ~iscell(props); props = {props}; end
    P = struct;
    for p = 1:numel(props)
        P.(props{p}) = [];
    end
    P.strain = [];
    P.rep = [];
    P.weight = [];
    
    tic;
    fprintf('Retrieving flat properties... \n')
    for d = 1:numel(obj.activeDirList)
        dataMP2 = obj.loadDataMP2(obj.activeDirList{d});
        Pd = dataMP2.flat_cellProps(props);
        Pd.strain = obj.strainMap.strain(dataMP2.strain);
        Pd.rep = obj.activeDirList{d};
        Pd.weight = ones(numel(Pd.(props{1})),1)*1./numel(Pd.(props{1})); 
        P = [P Pd]; %#ok<AGROW>
        
        fprintf('\t%d of %d done. Time elapsed = %5.5f sec. \n',d,numel(obj.activeDirList),toc)
    end
    fprintf('\n')
    P(1) = [];
    
    end
    
    
    function P = singleCC_cellProps(obj,props,groupByCell,catEndPts)
    %P = singleCC_cellProps(obj,props,groupByCell,catEndPts). groupByCell
    %is set to true by default, catEndPts is set to false.
    
    if ~iscell(props); props = {props}; end
    if nargin<3
        groupByCell = true;
        catEndPts = false;
    elseif nargin<4
        catEndPts = false;
    end
        
    P = struct;
    for p = 1:numel(props)
        P.(props{p}) = [];
    end
    P.strain = [];
    P.rep = [];
    P.weight = [];
    
    tic;
    fprintf('Retrieving single cell cycle properties...\n')
    for d = 1:numel(obj.activeDirList)
        dataMP2 = obj.loadDataMP2(obj.activeDirList{d});
        Pd = dataMP2.singleCC_cellProps(props,groupByCell,catEndPts);
        Pd.strain = obj.strainMap.strain(dataMP2.strain);
        Pd.rep = obj.activeDirList{d};
        Pd.weight = ones(numel(Pd.(props{1})),1)*1./numel(Pd.(props{1})); 
        P = [P Pd]; %#ok<AGROW>
        
        fprintf('\t%d of %d done. Time elapsed = %5.5f sec. \n',d,numel(obj.activeDirList),toc)
    end
    fprintf('\n')
    
    P(1) = [];
    
    end
    
    
    
    function Q = groupProps(obj,P,groupMthd,flatOrCC,normC)
    %Q = groupProps(obj,P,groupMethod,flatOrCC,normC)
    %Reshapes the structure P output by dataLabMP2.flat_cellProps, grouping
    %replicates by groupMthd = 'strain', 'loci', 'promoter', 'overnight',
    %or 'none' (the last simply flattens all data). Set normC to true to
    %normalize any concentrations or total fluorescence by their replicate
    %means. Set flatOrCC to 'flat' or 'CC' depending on whether P was
    %produced by flat_cellProps or singleCC_cellProps, this parameter only
    %required if normalization performed
    
    flds = fields(P);
            
    %Normalize any concentration related properties
    if nargin>4
        if normC~=0
            for f = 1:numel(flds)
                if (any(regexp(flds{f},'C')) || any(regexp(flds{f},'Itot'))) && ~any(regexp(flds{f},'growthRate')) && ~any(strfind(flds{f},'CM'))
                    switch flatOrCC
                        case 'flat'
                            for k = 1:numel(P)
                                P(k).(flds{f}) = P(k).(flds{f})./mean(P(k).(flds{f}));
                            end
                        case 'CC'
                            if normC==1
                                %Normalize using all data points
                                for k = 1:numel(P)
                                    Pk = P(k).(flds{f});
                                    nk = cellfun(@numel,Pk);
                                    Pk = cell2mat(Pk(:));
                                    Pk = Pk./mean(Pk);
                                    Pk = mat2cell(Pk,nk,1);
                                    P(k).(flds{f}) = Pk;
                                end
                            elseif normC==2
                                %Normalize using only the first data point
                                for k = 1:numel(P)
                                    Pk = P(k).(flds{f});
                                    meank = mean(cellfun(@(x) x(1),Pk));
                                    Pk = cellfun(@(x) x./meank,Pk,'uniformoutput',0);
                                    P(k).(flds{f}) = Pk;
                                end
                            elseif normC==3
                                %Normalize using each cell cycle's mean
                                for k = 1:numel(P)
                                    Pk = P(k).(flds{f});
                                    meank = mean(cellfun(@mean,Pk));
                                    Pk = cellfun(@(x) x./meank,Pk,'uniformoutput',0);
                                    P(k).(flds{f}) = Pk;
                                end
                            end                                
                    end
                end
            end
        end
    end
    
    %Default to no grouping if unspecified
    if nargin==2
        groupMthd = 'none';
    end
        
    strains = vertcat(P(:).strain);
    strainsUnique = unique(strains);

    Q = struct;
    switch groupMthd
        case 'strain'
            for s = 1:numel(strainsUnique)
                for p = 1:numel(flds)
                    if strcmp(flds{p},'rep')
                        Q(s).rep = {P(strains==strainsUnique(s)).rep}';
                    else
                        Q(s).(flds{p}) = vertcat(P(strains==strainsUnique(s)).(flds{p}));
                    end
                end
                Q(s).strain = Q(s).strain(1);
            end
        case 'loci'
            %Convert strains to loci
            loci = 0*strains;
            for s = 1:numel(strains)
                loci(s) = obj.strainMap.loci(obj.strainMap.inverse(strains(s)));
            end
            lociUnique = unique(loci);

            %Group
            for ll = 1:numel(lociUnique)
                for p = 1:numel(flds)
                    if strcmp(flds{p},'rep')
                        Q(ll).rep = {P(loci==lociUnique(ll)).rep}';
                    else
                        Q(ll).(flds{p}) = vertcat(P(loci==lociUnique(ll)).(flds{p}));
                    end
                end
                Q(ll).loci = lociUnique(ll);
            end

        case 'promoter'
            %Convert strains to loci
            pro = 0*strains;
            for s = 1:numel(strains)
                pro(s) = obj.strainMap.promoter(obj.strainMap.inverse(strains(s)));
            end
            proUnique = unique(pro);

            for pp = 1:numel(proUnique)
                for p = 1:numel(flds)
                    if strcmp(flds{p},'rep')
                        Q(pp).rep = {P(pro==proUnique(pp)).rep}';
                    else                        
                        Q(pp).(flds{p}) = vertcat(P(pro==proUnique(pp)).(flds{p}));
                    end
                end
                Q(pp).promoter = proUnique(pp);
            end

        case 'overnight'
            %Get replicates and determine the overnight
            overnight = 0*strains;
            for s = 1:size(P,2)
                overnight(s) = obj.strainMap.media(obj.getOvernight(P(s).rep));
            end

            overnightUnique = unique(overnight);
            for oo = 1:numel(overnightUnique)
                for p = 1:numel(flds)
                    if strcmp(flds{p},'rep')
                        Q(oo).rep = {P(overnight==overnightUnique(oo)).rep}';
                    else
                        Q(oo).(flds{p}) = vertcat(P(overnight==overnightUnique(oo)).(flds{p}));
                    end
                end
                Q(oo).overnight = overnightUnique(oo);
            end

        case 'none'
            %Simply flatten the entire structure
            for p = 1:numel(flds)
                if strcmp(flds{p},'rep')
                    Q.rep = {P(:).rep}';
                else
                    Q.(flds{p}) = vertcat(P(:).(flds{p}));
                end
            end

    end
    
    end
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% NUMERIC CODING FOR STRAINS, LOCI, PROMOTERS, MEDIA

    function C = initStrainMap(obj)
        C = struct;
        C.strain = containers.Map;
        C.loci = containers.Map;
        C.promoter = containers.Map;
        
        strainN = 1:numel(obj.strainList);
        lociN = strainN(end) + (1:numel(obj.lociList));
        promoterN = lociN(end) + (1:numel(obj.promoterList));
        
        for s = 1:numel(obj.strainList)
            C.strain(obj.strainList{s}) = strainN(s);
            C.loci(obj.strainList{s}) = lociN(strcmp(obj.lociList,obj.strainParse(obj.strainList{s},'loci')));
            C.promoter(obj.strainList{s}) = promoterN(strcmp(obj.promoterList,obj.strainParse(obj.strainList{s},'promoter')));
        end
        
        mediaN = promoterN(end) + (1:numel(obj.mediaList)+1);
        C.media = containers.Map([obj.mediaList 'LB'], mediaN);
    
        C.inverse = containers.Map(1:mediaN(end),[obj.strainList' obj.lociList obj.promoterList obj.mediaList 'LB']);
    end
    
    function str = strainParse(obj,strain,mthd)
        switch mthd
            case 'loci'
                for k = 1:numel(obj.lociList)
                    if any(strfind(strain,obj.lociList{k}))
                        str = obj.lociList{k};
                        return
                    end
                end
            case 'promoter'
                str = cell(numel(obj.promoterList),1);
                for k = 1:numel(obj.promoterList)
                    if any(strfind(strain,obj.promoterList{k}))
                        str{k} = obj.promoterList{k};
                    end
                end
                nchar = cellfun(@numel,str);
                str = str(nchar==max(nchar));
                str = str{1};
        end
    end
    
    function str = inverseCode(obj,v)
    %str = inverseCode(obj,v). Takes a numeric vector and returns a cell of
    %character strings with the inverse code
    str = cell(numel(v),1);
    for n = 1:numel(v)
        str{n} = obj.strainMap.inverse(v(n));
    end
    
    end
    
    function R = initRGB(obj)
    omitStrains = {'P5left3','gyrAori3'};
    strainListShort = obj.strainList(~ismember(obj.strainList,omitStrains));
    strainListShort = strainListShort(:);
    keys = [strainListShort' obj.lociList obj.promoterList];
    rgb = [sqrt(0.5*hsv(numel(strainListShort))+0.1); 0.85*copper(numel(obj.lociList)); 0.75*jet(numel(obj.promoterList))];
    rgb = mat2cell(rgb,ones(size(rgb,1),1),3);        
    R = containers.Map(keys,rgb);
    
    for k = 1:numel(omitStrains)
        R(omitStrains{k}) = [0 0 0];
    end
    R('M9S') = [0.2 0.8 1];
    R('M9F') = [1 0.8 0.2];
    R('LB') = [1 0 0.5];
    
    end
    
    function rgb = getRGB(obj,str)
    str = char(str);
    if isKey(obj.RGB,str)
        rgb = obj.RGB(str);
    elseif any(strfind(str,'deadEnd'))
        baseProp = regexp(str,'(%\w+)deadEnd','tokens');
        if ~isempty(baseProp)
            baseProp = baseProp{1}{1};
            if isKey(obj.RGB,baseProp)
                rgb = sqrt(obj.RGB(baseProp));
            else
                rgb = [0 0 0];
            end
        else
            rgb = [0 0 0];
        end
    else
        rgb = [0 0 0];                
    end
    
    end
    
    
    
    function hf = rgbLegends(obj,str,textColor)
    %Outputs a figure mapping color to the strains, loci, or promoters in
    %string
    if ~iscell(str), str = {str}; end
    if nargin<2, textColor = 'k'; end
    dy = 1/numel(str);
    
    hf = figure;
    ha = axes('parent',hf,'position',[0 0 1 1]);
    hold(ha,'on')
    for k = 1:numel(str)
        xx = [0 0 1 1 0];
        yy = 1+dy*[-k -(k-1) -(k-1) -k -k];
        fill(xx,yy,obj.getRGB(str{k}),'edgecolor','w','parent',ha,'linewidth',2)
        text(0.5,1-dy*(k-0.5),str{k},'color',textColor,'horizontalalignment','center','fontsize',24)
    end
    
        
    end
    
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% MISC FUNCTIONS
    
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
    
    
    
    
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    methods (Static = true)
    
    function overnight = getOvernight(baseDir)
        %overnight = GETOVERNIGHT(obj,baseDir)
        yy = baseDir(1:4);
        if str2double(yy) < 2016
            overnight = 'LB';
        else
            overnight = 'M9S';
        end
    end
    
    function str = reshapeStr(str0,nCharsWide)
        str = cell(ceil(numel(str0)/nCharsWide),1);
        strLeft = str0;
        
        c = 1;
        while any(strLeft)
            strc = strLeft(1:min(nCharsWide-1,numel(strLeft)));
            whiteSpace = strfind(strc,' ');
            if ~any(whiteSpace) || all(whiteSpace==1)
                whiteSpace = numel(strc);
            end
            str{c} = strc(1:whiteSpace(end));
            strLeft = strLeft(whiteSpace(end)+1:end);
            c = c + 1;
        end
        
    end
    
    function P2 = flat_cellPropsFilt(P,filterNames,filterRange)
    %P2 = flat_cellPropsFilt(P,filterNames,filterRange)
    %Using the structure output by dataLabMP2.flat_cellProps (or
    %dataMP2.flat_cellProps), filters out data according to filterNames and
    %and filterRange. 
    if ~iscell(filterNames); filterNames = {filterNames}; end
    if ~iscell(filterRange); filterRange = {filterRange}; end
    
    flds = fields(P);
    useMe = false(numel(filterNames),1);
    for f = 1:numel(filterNames)
        if any(strcmp(flds,filterNames{f}))
            useMe(f) = true;
        else
            warning('%s is not an accessible field. Removed.',filterNames{f})
        end
    end
    
    if sum(useMe)==0
        warning('No filtering could be performed.')
        P2 = P;
        return
    end
    
    filterNames = filterNames(useMe);
    filterRange = filterRange(useMe);
        
    P2 = P;
    flds = fields(P);
    for p = 1:numel(P)
        nPts = numel(P(p).(filterNames{1}));
        Z = true(nPts,1);
        for f = 1:numel(filterNames)
            Z = Z & (P(p).(filterNames{f}) >= filterRange{f}(1)) & (P(p).(filterNames{f}) < filterRange{f}(2)); 
        end
        
        for f = 1:numel(flds)
            canFilter = numel(P(p).(flds{f}))==nPts;
            if canFilter
                P2(p).(flds{f}) = P(p).(flds{f})(Z);
            end
        end
    end
    
    end

    
    
    end
end