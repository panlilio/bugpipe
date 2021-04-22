function evalGUIparent(D,fov,startFrame,L)

if ispc
    hf = figure('position',[1537,409,955,832]);
else
    hf = figure('position',[264 42 1020 953]);
end

if nargin==3
    L = imListMP(D.dataDir);
end

frames = L.imListFOVT(L.imListFOVT(:,1)==fov,2);

currentFrame = startFrame;

initGUI;

function initGUI
    if ~exist(fullfile(D.mat2Dir,sprintf('fov%.2d_t%.4d.mat',fov,currentFrame)),'file')
        disp(D.mat2Dir)
        fprintf('fov%.2d_t%.4d.mat does not exist. Try a different frame/fov.\n',fov,currentFrame)
        return
    end
    
    clf(hf)
    
    C0 = load(fullfile(D.mat2Dir,sprintf('fov%.2d_t%.4d.mat',fov,currentFrame)));
    C1 = load(fullfile(D.mat2Dir,sprintf('fov%.2d_t%.4d.mat',fov,currentFrame+1)));
    if nargin==3
        L = imListMP(D.dataDir);
    end
    
    I0 = imread(L.getfilenamesFOVT(fov,currentFrame,1));
    I0 = double(I0);
    
    I1 = imread(L.getfilenamesFOVT(fov,currentFrame+1,1));
    I1 = double(I1);
    
    t0 = currentFrame;
    t = currentFrame+1;
    evalGUI0(C0,C1,I0,I1,t0,t,hf);
    
    uicontrol('parent',hf,'style','pushbutton','string','<html>PREVIOUS<br>FRAME','units','normalized','position',[0.74 0.35 0.1 0.08],'callback',@prevFrameCallback)
    uicontrol('parent',hf,'style','pushbutton','string','<html>NEXT<br>FRAME','units','normalized','position',[0.85 0.35 0.1 0.08],'callback',@nextFrameCallback)
    uicontrol('parent',hf,'style','edit','string',num2str(currentFrame),'units','normalized','position',[0.85 0.29 0.04 0.04],'callback',@pickFrameCallback)
    uicontrol('parent',hf,'style','text','string','JUMP TO FRAME','units','normalized','horizontalalignment','right','position',[0.74 0.28 0.1 0.04])
    
end

function prevFrameCallback(source,eventData)
    if ismember(currentFrame-1,frames)
        currentFrame = currentFrame - 1;
        initGUI
    else
        fprintf('Frame = %d cannot be loaded. \n',currentFrame-1)
    end
end

function nextFrameCallback(source,eventData)
    if ismember(currentFrame + 1,frames)
        currentFrame = currentFrame + 1;
        initGUI
    else
        fprintf('Frame = %d cannot be loaded. \n',currentFrame+1)
    end
end

function pickFrameCallback(source,eventData)
    str = get(source,'string');
    val = str2double(str);
    if ismember(val,frames)
        currentFrame = val;
        initGUI
    else
        fprintf('Frame = %d cannot be loaded. \n',val)
    end
end

end