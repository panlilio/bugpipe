function dirList = getDirsMP_SSD(baseDir)
%DIRLIST = GETDIRSMP_SSD(BASEDIR). 
%Retrieves a list of directories containing data relevant to the experiment
%of interest. BASEDIR only needs the experiment [date], not a full path.
%Builds directories to write MATs to SSD (C:).
%OUTPUT FIELDS:
% .baseDir  experiment date
% .dataDir  pngs
% .matDir   mats with CELLS and CHANNELS from segmentation
% .mat2Dir  mats with CELLS and CHANNELS from tracking
% .lineDir  mats with lineage info
% .textDir  txts, flat
% .textLineDir  txts of lineages, subdirectory of textDir

if nargin==0 
    baseDir = uigetdir('M:\Data_png');
    s = strfind(baseDir,filesep);
    baseDir = baseDir(s(end)+1:end);
end

dataInfo = lookupDataMP(baseDir);

dirList = {};
dirList.baseDir = baseDir;
dirList.dataDir = fullfile('M:\Data_png',dataInfo.STRAIN,baseDir);
if exist([dirList.dataDir '_sorted'],'dir')
    dirList.dataDir = [dirList.dataDir '_sorted'];
end
dirList.matDir = fullfile('C:\Users\Mia\Desktop\Data_temp\Data_mat',dataInfo.STRAIN,baseDir);
dirList.mat2Dir = fullfile('C:\Users\Mia\Desktop\Data_temp\Data_mat2',dataInfo.STRAIN,baseDir);
dirList.lineDir = fullfile('C:\Users\Mia\Desktop\Data_temp\Data_line',dataInfo.STRAIN,baseDir);
dirList.textDir = fullfile('C:\Users\Mia\Desktop\Data_temp\Data_text',dataInfo.STRAIN,baseDir);
dirList.textLineDir = fullfile(dirList.textDir,'lineages');

if ispc
    dirList.flatDir = fullfile('C:\Users\Mia\Desktop\Data_temp\Data_flat',dataInfo.STRAIN,baseDir);
else
    dirList.flatDir = fullfile('/Users/miapanlilio/Dropbox/MP/dataLabMP/data/Data_flat',dataInfo.STRAIN,baseDir);
end

end