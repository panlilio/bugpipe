classdef imListMP < handle
    properties
        pngDir
        imList
        blankList
        imListFOVT
        blankListFOVT
        timeCreated
        fovRange
        tRange
    end

    methods
        function obj = imListMP(pngDir)            
            obj.pngDir = pngDir;
            
            pngList = ls(pngDir);
            pngList = pngList(3:end,:);
            pngList = cellstr(pngList);
            
            notHeader = cellfun(@(x) ~any(strfind(x,'header.mat')),pngList);
            isBlank = cellfun(@(x) any(strfind(x,'BLANK')),pngList);
            isBlank = isBlank & notHeader;
            
            isFluor = cellfun(@(x) any(strfind(x,'ch01')),pngList);
                        
            imList = pngList(~isBlank & notHeader & isFluor);
            blankList = pngList(isBlank & notHeader & isFluor);
            
            fovT = regexp(imList,'(\d*)_t(\d*)','tokens');
            [fov,t] = cellfun(@(x) x{:}{:},fovT,'uniformoutput',0);
            fovT = [str2double(fov) str2double(t)];
            
            obj.imList = imList;
            obj.imListFOVT = fovT;
            
            fovT = regexp(blankList,'(\d*)_t(\d*)','tokens');
            [fov,t] = cellfun(@(x) x{:}{:},fovT,'uniformoutput',0);
            fovT = [str2double(fov) str2double(t)];
            
            obj.blankList = blankList;
            obj.blankListFOVT = fovT;
            
            actualTime = regexp(imList,'[.](\d{2}\w*\d{4})_(\d{2})[.](\d{2})[.](\d{2}).png','tokens');
            [dmy,hh,mm,ss] = cellfun(@(x) x{:}{:},actualTime,'uniformoutput',0);   
            hh = str2double(hh);
            mm = str2double(mm);
            ss = str2double(ss);
            
            TT = datenum(dmy);
            
            TT = TT*24*60 + hh*60 + mm + ss/60;
            TT = TT - TT(1);
            
            if strcmp(obj.pngDir(end-7:end),'20160325')
                %DAYLIGHT SAVINGS TIME: modify to account for hour jump
                %at 1 am.
                TT = obj.dstJump(TT);                
            end
            
            obj.timeCreated = TT;    
            
            obj.fovRange = unique(obj.imListFOVT(:,1))';
            obj.tRange = unique(obj.imListFOVT(:,2))';
            
            blankList = cell(size(imList));
            for k = 1:numel(imList)
                [~,blankname] = obj.getfilenames(k,0);
                blankList{k} = blankname;
            end
            obj.blankList = blankList;
            obj.blankListFOVT = obj.imListFOVT;
        end
        
        function [imname,blankname] = getfilenames(obj,n,includeParent)    
            imname = obj.imList{n};
            if includeParent
                imname = fullfile(obj.pngDir,imname);
            end
            fov = obj.imListFOVT(n,1);
            t = obj.imListFOVT(n,2);
            [~,blankname] = obj.getfilenamesFOVT(fov,t,includeParent);
        end
        
        function [imname,blankname] = getfilenamesFOVT(obj,fov,t,includeParent)
            imID = obj.imListFOVT(:,1)==fov & obj.imListFOVT(:,2)==t;
            if sum(imID)==0
                imname = [];
            elseif sum(imID)>1
                imname = obj.imList(imID);
                imname = imname{end};
                if includeParent
                    imname = fullfile(obj.pngDir,imname);
                end
            else
                imname = obj.imList{imID};
                if includeParent
                    imname = fullfile(obj.pngDir,imname);
                end
            end
            
            blankID = obj.blankListFOVT(:,1)==fov & obj.blankListFOVT(:,2)==t;
            if sum(blankID)==0
                blankname = [];
            elseif sum(blankID)>1
                blankname = obj.blankList(blankID);
                blankname = blankname{end};
            else
                blankname = obj.blankList{blankID};
            end
            
            if ~isempty(blankname)
                if includeParent
                    blankname = fullfile(obj.pngDir,blankname);
                end
            end
        end
        
        function tAcq = tAcq(obj,fov,t)
            tID = obj.imListFOVT(:,1)==fov & obj.imListFOVT(:,2)==t;
            tAcq = obj.timeCreated(tID);
        end
    end
    
    
    methods (Hidden = true, Static = true)
        function TT = dstJump(TT)
        %Daylight savings time: adjust time acquired after 1am (1 h jump)      
            TTsorted = sort(TT);
            dTTsorted = diff(TTsorted);
            jumped = find(dTTsorted==max(dTTsorted));
            timeJumped = TTsorted(jumped+1);
            TT(TT>=timeJumped) = TT(TT>=timeJumped) - 60;
        end
    end
end