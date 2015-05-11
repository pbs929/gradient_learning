function [M, MClean, MNoise] = loadMat_noise(param,type)
% Loads multiple speech files, computes peripheral response and delays, 
% and concatenates into a single data matrix
% type = 'train' or 'test' specifies which input folder to use
% This version loads multiple noise mixes as well.

switch type
    case 'train'
        dataDir = param.trainDir;
        whichFiles = param.trainFiles;
    case 'test'
        dataDir = param.testDir;
        whichFiles = param.testFiles;
    otherwise 
        error('Invalid data type option - use ''train'' or ''test''');
end

% Clean data all together
M = dataLoad(dataDir, '*clean.mat', whichFiles, param);
if nargout==1, return; end

% Clean data separated by noise mix
MClean = cell(1,length(param.noiseTypes));
for iNse = 1:length(param.noiseTypes)
    nse = param.noiseTypes(iNse);
    MClean{iNse} = dataLoad(dataDir, ['*clean',num2str(nse),'.mat'], whichFiles, param);
end
if nargout==2, return; end

% Noisy data
MNoise = cell(1,length(param.noiseTypes));
for iNse = 1:length(param.noiseTypes)
    nse = param.noiseTypes(iNse);
    MNoise{iNse} = dataLoad(dataDir, ['*n',num2str(nse),'_snr',num2str(param.noiseLevs),'.mat'], whichFiles, param);
end

% Noise mixes 

function M = dataLoad(dataDir, name, whichFiles, param)

% get list of speech files
files = rdir([dataDir, name]);
if ~isempty(whichFiles), files = files(whichFiles); end

% loop through files and append to data matrix 
M = [];
for fle = 1:length(files)
    [~,~,~,data]=load_data(files(fle).name,dataDir,param.periph,...
        param.normType,param.int,param.delay);
    data = data(:,param.delay:end);
    disp(files(fle).name)
    
    M=[M,data];
end

