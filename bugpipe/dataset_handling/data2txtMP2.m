function data2txtMP2(data_or_lab,varargin)
%data2txtMP2(data_or_lab,propsA,propsB,destdir)
%Takes input data (dataMP2) or lab (dataLabMP2) and writes properties to
%file in destdir. 
%
%propsA     cell array of strings containing all properties to be
%           written in addition to 'cellID', 'parentID', 'line', 'fov'
%
%propsB     the same as propsA but for flat files, i.e. start and end
%           points
%
%destdir    destination directory to which data will be written. If no 
%           input is given, the current working directory is used

%%%%%%%%%%%%%%%%%%%%%%% PARSE INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = inputParser;
addRequired(p,'data_or_lab',@(x) isa(x,'dataLabMP2') || isa(x,'dataMP2'))
addOptional(p,'propsA',{},@iscell)
addOptional(p,'propsB',{},@iscell)
addOptional(p,'destdir',cd,@ischar)
parse(p,data_or_lab,varargin{:});

%%%%%%%%%%%%%%%%%%%%%%% WRITE TO FILE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isa(data_or_lab,'dataMP2')
    data_or_lab.writedata(p.Results.propsA,p.Results.propsB,destdir)
else
    %Write an additional readme for a lab with more than one replicate
    readme_name = fullfile(p.Results.destdir,'data_readme.txt');
    fid = fopen(readme_name,'w');
    fprintf(fid,'DATE: %s\n\n',date);
    fprintf(fid,'FILTERS USED \n');
    fiu = fields(data_or_lab.dataFilter);
    for ff = 1:numel(fiu)
        if isempty(data_or_lab.dataFilter.(fiu{ff}))
            continue
        elseif ischar(data_or_lab.dataFilter.(fiu{ff}))
            fprintf(fid,'\t%s: %s\n',fiu{ff},data_or_lab.dataFilter.(fiu{ff}));
        elseif numel(data_or_lab.dataFilter.(fiu{ff}))==1 
            fprintf(fid,'\t%s: %d\n',fiu{ff},data_or_lab.dataFilter.(fiu{ff}));
        else
            fprintf(fid,'\t%s: [%f,%f]\n',fiu{ff},data_or_lab.dataFilter.(fiu{ff}));
        end
    end
    fprintf(fid,'\nSTRAIN \tREPLICATE \tORIENTATION \tSWITCH \tMEDIA \tOVERNIGHT \tNCELLS \tCOMMENTS \n\n');
    fclose(fid);
    
    for k = 1:numel(data_or_lab.activeDirList)
        ak = data_or_lab.loadDataMP2(data_or_lab.activeDirList{k});
        [~,nk] = ak.writedata(p.Results.propsA,p.Results.propsB,p.Results.destdir);
        
        shiftT = num2str(ak.switchTime,'%d,');
        shiftT = shiftT(1:end-1);
        fid = fopen(readme_name,'a');
        fprintf(fid,'%s \t%s \t%s \t%s \t%s \t%s \t%d \t%s \n',...
            ak.strain, ak.dirList.baseDir, ak.orientation, shiftT, ak.media, ak.getOvernight, nk, ak.comments);
        fclose(fid);
        
        fprintf('%d of %d replicates written to file.\n',k,numel(data_or_lab.activeDirList));
    end
    
    
end

end