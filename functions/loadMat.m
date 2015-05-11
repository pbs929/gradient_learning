function A = loadMat(param)
% Loads multiple speech files, computes peripheral response and delays, 
% and concatenates into a single data matrix

% get list of speech files
files = rdir([param.sampleDir, '*.mat']);
if ~isempty(param.whichFiles), files = files(param.whichFiles); end

% loop through files and append to data matrix 
A = [];
for fle = 1:length(files)
    [~,~,~,data]=load_data(files(fle).name,param.sampleDir,param.periph,...
        param.normType,param.int,param.delay);
    data = data(:,param.delay:end);
    disp(files(fle).name)
    
    A=[A,data];
end