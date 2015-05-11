function [M, MClean, MNoise] = loadMat_gauss(nRows, nCols, param)
% Loads random gaussian distributed data in the formats required by the
% GradDesc script.  

% Clean data all together
M = randn(nRows,nCols);
if nargout==1, return; end

% Clean data separated by noise mix
MClean = cell(1,length(param.noiseTypes));
for iNse = 1:length(param.noiseTypes)
    MClean{iNse} = randn(nRows,nCols);
end
if nargout==2, return; end

% Noisy data
MNoise = cell(1,length(param.noiseTypes));
for iNse = 1:length(param.noiseTypes)
    MNoise{iNse} = randn(nRows,nCols);
end
