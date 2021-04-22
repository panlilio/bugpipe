function dataInfo = lookupDataMP(directoryName)

if nargin==0 
    directoryName = uigetdir('M:\Data_png');
end

s = strfind(directoryName,filesep);
if any(s)
    directoryName = directoryName(s(end)+1:end);
end

if ispc
    fid = fopen('C:\Users\Mia\Desktop\dataBook.txt','r');
else
    fid = fopen('/Users/miapanlilio/Desktop/dataLabMP_standalone/dataBook.txt','r');
end

nCols = 9;
strFormat = repmat('%s',1,nCols);
strFormat = [strFormat '\n'];
headers = textscan(fid,strFormat,'delimiter','\t');
strFormat = strFormat(1:end-2);
data0 = textscan(fid,strFormat,'emptyvalue',NaN,'delimiter','\t','whitespace','');

f = @(x) x{:};
headers = cellfun(f,headers,'uniformoutput',false);

s = strcmp(headers,'FOLDER');
s = strcmp(strtrim(data0{s}),directoryName);

data = [];
for c = 1:nCols
    dataC = strtrim(data0{c}(s));
    dataC = dataC{:};
    switch dataC
        case {'Y','y','1'}
            dataC = true;
        case {'N','n','0'}
            dataC = false;
        case {'L','l'}
            dataC = 'left';
        case {'R','r'}
            dataC = 'right';
        case {'NA','Na','na'}
            dataC = '0';
    end
    dataC = {dataC};
    data = [data dataC];
end

dataInfo = cell2struct(data,headers,2);
dataInfo.SWITCH = str2num(dataInfo.SWITCH); %#ok<ST2NM>
dataInfo.FRAMELIM = str2num(dataInfo.FRAMELIM); %#ok<ST2NM>

fclose all;

end