function hf = evalGUI0(C0,C1,I0,I1,t0,t,hf)
%Assess variables to initialize figure
C = {};
I = {};
I{1} = I0;
I{2} = I1;
if ~isempty(C0)
    C{1} = C0;
    C{2} = C1;
    cLab = unique([C0.CHANNELS.label C1.CHANNELS.label]);
    nc = numel(cLab);
    rgbChannels = containers.Map(cLab,mat2cell(hsv(nc),ones(1,nc),3));
    cellLab = unique([C0.CELLS.label; C1.CELLS.label]);
    nCells = numel(cellLab);
    rgbCells = containers.Map(cellLab,mat2cell(hsv(nCells),ones(1,nCells),3));
else
    C{2} = C1;
    cLab = unique(C1.CHANNELS.label);
    nc = numel(cLab);
    rgbChannels = containers.Map(cLab,mat2cell(hsv(nc),ones(1,nc),3));
    cellLab = unique(C1.CELLS.label);
    nCells = numel(cellLab);
    rgbCells = containers.Map(cellLab,mat2cell(hsv(nCells),ones(1,nCells),3));
end

%Initialize figure
% hf = figure('position',[260 6 1020 990]);

%Panel for fluorescent images
hpImages = uipanel(hf,'units','normalized','position',[0.005 0.505 0.99 0.49]);
haImages = zeros(1,2);
haImages(1) = axes('parent',hpImages,'units','normalized','position',[0.005 0.005 0.49 0.99]);
haImages(2) = axes('parent',hpImages,'units','normalized','position',[0.505 0.005 0.49 0.99]);

%Panel for specific channel
hpChannels = uipanel(hf,'units','normalized','position',[0.005 0.25 0.49 0.25],'title','SELECT CHANNEL');
haChannels = zeros(1,3);
haChannels(1) = axes('parent',hpChannels,'position',[0.005 0.43 0.99 0.42],'xtick','','ytick','');
haChannels(2) = axes('parent',hpChannels,'position',[0.005 0.005 0.99 0.42],'xtick','','ytick','');
selectChannel = uicontrol('parent',hpChannels,'style','popupmenu','string',cellstr(num2str(cLab(:))),...
    'units','normalized','position',[0.3 0.8 0.4 0.19],'callback',@selectChannelCallback);

%Panel with cell properties from specific channel
hpCellProps = uipanel(hf,'units','normalized','position',[0.005 0.005 0.99 0.24],'title','CELL PROPS');
htCellProps = uitable('parent',hpCellProps,'units','normalized','position',[0.005 0.005 0.99 0.99]);

%Panel with frame-to-frame stats
hpStats = uipanel(hf,'units','normalized','position',[0.505 0.25 0.49 0.25],'title','FRAME-TO-FRAME STATS');
hDivText = uicontrol('parent',hpStats,'style','text','string','DIVISION EVENTS','units','normalized',...
    'position',[0.1 0.8 0.2 0.1]);
hChannelsLostText = uicontrol('parent',hpStats,'style','text','string','CHANNELS LOST','units','normalized',...
    'position',[0.1 0.65 0.2 0.1]);
hChannelsFoundText = uicontrol('parent',hpStats,'style','text','string','NEW CHANNELS','units','normalized',...
    'position',[0.1 0.55 0.2 0.1]);
hCellsLostText = uicontrol('parent',hpStats,'style','text','string','CELLS LOST','units','normalized',...
    'position',[0.1 0.4 0.2 0.1]);
hCellsFoundText = uicontrol('parent',hpStats,'style','text','string','NEW CELLS','units','normalized',...
    'position',[0.1 0.3 0.2 0.1]);


n = divCount(C0,C1);
hDivTextN = uicontrol('parent',hpStats,'style','text','string',num2str(n),'units','normalized',...
    'position',[0.3 0.8 0.1 0.1]);

n = channelCount(C0,C1);
hChannelsLostTextN = uicontrol('parent',hpStats,'style','text','string',num2str(n(1)),'units','normalized',...
    'position',[0.3 0.65 0.1 0.1]);
hChannelsFoundTextN = uicontrol('parent',hpStats,'style','text','string',num2str(n(2)),'units','normalized',...
    'position',[0.3 0.55 0.1 0.1]);

n = cellCount(C0,C1);
hCellsLostTextN = uicontrol('parent',hpStats,'style','text','string',num2str(n(1)),'units','normalized',...
    'position',[0.3 0.4 0.1 0.1]);
hCellsFoundTextN = uicontrol('parent',hpStats,'style','text','string',num2str(n(2)),'units','normalized',...
    'position',[0.3 0.3 0.1 0.1]);


%DISPLAY IMAGES AND LABELLED CELLS
%Fluorescent images
imshow(mat2gray(I0),'parent',haImages(1))
imshow(mat2gray(I1),'parent',haImages(2))

hold(haImages(1),'on')
hold(haImages(2),'on')

text(5,12,sprintf('frame = %g',t0),'parent',haImages(1),'color','w','fontsize',12)
text(5,30,sprintf('t = %6.2f min',C{1}.CELLS.tAcq),'parent',haImages(1),'color','w','fontsize',10)
text(5,48,sprintf('N = %d',numel(C{1}.CELLS.lineID)),'parent',haImages(1),'color','w','fontsize',8)
text(5,12,sprintf('frame = %g',t),'parent',haImages(2),'color','w','fontsize',12)
text(5,30,sprintf('t = %6.2f min',C{2}.CELLS.tAcq),'parent',haImages(2),'color','w','fontsize',10)
text(5,48,sprintf('N = %d',numel(C{2}.CELLS.lineID)),'parent',haImages(2),'color','w','fontsize',8)


for k = 1:numel(C)
    if isempty(C{k}), continue, end
    
    %Plot and label channels
    for j = 1:numel(C{k}.CHANNELS.threshold)
        bb = C{k}.CHANNELS.boundingBox(j,:);
        xx = [bb(3) bb(3) bb(4) bb(4) bb(3)];
        yy = [bb(1) bb(2) bb(2) bb(1) bb(1)];
        plot(xx,yy,'color',rgbChannels(C{k}.CHANNELS.label(j)),'parent',haImages(k))
        text(xx(1)-20,yy(1),num2str(C{k}.CHANNELS.label(j)),'color',rgbChannels(C{k}.CHANNELS.label(j)),'parent',haImages(k))
        plotCells(C{k},C{k}.CHANNELS.label(j),[1 1 1 1],haImages(k),0)
    end
    
end

%% CALLBACK FUNCTIONS
function selectChannelCallback(source,eventdata)
    channelstr = get(source,'string');
    channel = get(source,'value');
    channel = str2double(channelstr(channel));
    
    if ~isempty(C0)
        if any(C0.CHANNELS.label==channel)
            cla(haChannels(1))
            chID = C0.CHANNELS.label==channel;
            bbch = C0.CHANNELS.boundingBox(chID,:);
            Iselect = I0(bbch(1):bbch(2),bbch(3):bbch(4));
            imshow(mat2gray(Iselect),'parent',haChannels(1))
            plotCells(C0,channel,bbch,haChannels(1),1)
            axis(haChannels(1),'image')
        end
    end
    
    if any(C1.CHANNELS.label==channel)
        cla(haChannels(2))
        chID = C1.CHANNELS.label==channel;
        bbch = C1.CHANNELS.boundingBox(chID,:);
        Iselect = I1(bbch(1):bbch(2),bbch(3):bbch(4));
        imshow(mat2gray(Iselect),'parent',haChannels(2))
        plotCells(C1,channel,bbch,haChannels(2),1)
        axis(haChannels(2),'image')
    end
    
    T = getCellStats(C,channel,I);
    set(htCellProps,'Data',T.data,'RowName',T.rows,'ColumnName',T.cols)
end

%% MISC FUNCTIONS
function plotCells(C,chLabel,bb,ha,doLabel)
    hold(ha,'on')
    cID = find(C.CELLS.channel==chLabel);
    cID = reshape(cID,1,numel(cID));
    n_ = 0.75;
    for cc = cID
        plot(C.CELLS.boundary{cc}(:,1)-bb(3)+1,C.CELLS.boundary{cc}(:,2)-bb(1)+1,'color',rgbCells(C.CELLS.label{cc}),'parent',ha)        
        if doLabel
            yf = (n_+0.2*rand(1))*(bb(2)-C.CELLS.centroid(cc,2)) + C.CELLS.centroid(cc,2);
            yf = [C.CELLS.centroid(cc,2) yf] - bb(1) + 1;
            xf = C.CELLS.centroid(cc,1);
            xf = [xf xf] - bb(3) + 1;
            plot(xf,yf,'color',rgbCells(C.CELLS.label{cc}),'parent',ha)
            %For full string label
%             cellstr = C.CELLS.label{cc};
%             s = strfind(cellstr,'_');
%             s = s(2);
%             cellstr = cellstr(s+1:end);
%             s = strfind(cellstr,'_');
%             cellstr(s) = '-';
%             text(xf(1),yf(2),cellstr,'color',rgbCells(C.CELLS.label{cc}),'parent',ha,'horizontalalignment','center','fontsize',8)
            
            %Simple lineID label
            text(xf(1),yf(2),num2str(C.CELLS.lineID(cc)),'color',rgbCells(C.CELLS.label{cc}),'parent',ha,'horizontalalignment','center','fontsize',8)
            n_ = -n_;
        end
    end
    
    if doLabel
        xxch = 1;
        yych = 3;
        cellstr = C.CELLS.label{cc};
        s = strfind(cellstr,'_');
        cellstr(s) = '-';
        text(xxch,yych,cellstr(1:s(2)-1),'color',rgbChannels(chLabel),'parent',ha,'fontsize',12,'fontweight','bold')
    end
end

function T = getCellStats(C,chLabel,I)
    Tdata = [];
    Trows = {};
    ac = 1;
    for CC = 1:numel(C)
        if ~isempty(C{CC})
            cID = find(C{CC}.CELLS.channel==chLabel);
            cID = reshape(cID,1,numel(cID));
            for cc = cID
                A0 = getCellPropsMP(C{CC}.CELLS.boundary{cc},I{CC});
                Tdata = [Tdata; A0.area A0.majorAxis A0.minorAxis A0.regionpropsCentroid A0.intensityMeanSD(1)/A0.area A0.intensityMeanSD(1) A0.intensityRange];
                if CC==1
                    Trows{ac} = sprintf('t=%g | %s | lineID = %d',t0,C{CC}.CELLS.label{cc}(9:end),C{CC}.CELLS.lineID(cc));
                else
                    Trows{ac} = sprintf('t=%g | %s | lineID = %d',t,C{CC}.CELLS.label{cc}(9:end),C{CC}.CELLS.lineID(cc));
                end
                ac = ac + 1;
            end
        end
    end
    Tcols = {'A' 'L' 'w' 'x' 'y' 'I/A' 'I' 'Imin' 'Imax' 'QC'};
    T.data = Tdata;
    T.rows = Trows;
    T.cols = Tcols;
end

function n = divCount(C0,C1)
    if isempty(C0)
        n = NaN;
    else
        n = 0;
        C0lab = C0.CELLS.label;
        C1lab = C1.CELLS.label;
        for cc = 1:numel(C0lab)
            s1 = strfind(C1lab,C0lab{cc});
            s1 = cellfun(@any,s1);
            s1 = s1(:);
            s2 = strcmp(C1lab,C0lab{cc});
            s3 = not(s2) & s1;
            if any(s3)
                n = n + 1;
            end
        end
    end
end


function n = channelCount(C0,C1)
    if isempty(C0)
        n = [NaN NaN];
    else
        n = [0 0];
        n(1) = numel(C0.CHANNELS.label) - sum(ismember(C0.CHANNELS.label,C1.CHANNELS.label));
        n(2) = numel(C1.CHANNELS.label) - sum(ismember(C1.CHANNELS.label,C0.CHANNELS.label));
    end
end

function n = cellCount(C0,C1)
    if isempty(C0)
        n = [NaN NaN];
    else
        n = [numel(C0.CELLS.label) 0];
        C0lab = C0.CELLS.label;
        C1lab = C1.CELLS.label;
        for cc = 1:numel(C0lab)
            s1 = strfind(C1lab,C0lab{cc});
            s1 = cellfun(@any,s1);
            s1 = s1(:);
            
            if any(s1)
                n(1) = n(1) - 1;
                n(2) = n(2) + sum(s1);
            end
        end
        n(2) = numel(C1lab) - n(2);
    end
end


end