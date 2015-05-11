function param = params_GradDesc()
% This function sets parameters for the Gradient descent experiment.  To 
% be used with the script GradDesc_AURORA.  PBS 2-14-14

% Data file selection
%param.trainDir = '../data/samples_aurora/data_aurora_train/'; 
%param.trainFiles = [1:5,56:60]; % male files are 56:109
param.trainDir = '../data/samples_aurora/data_aurora_test/'; 
param.trainFiles = [1:5,53:57]; % male files are 53:103 
param.testDir = '../data/samples_aurora/data_aurora_test/'; 
param.testFiles = [1:5,53:57]; % male files are 53:103

% Noise parameters for testing and training data
param.noiseTypes = 1:4;      % 4 types total
param.noiseLevs = 5;         % -5:5:20 dB available

% Peripheral model options
param.periph = 'gamma';      % spec-spectrogram gamma-gammatone meddis-meddis model NSL-NSL spectrogram
param.int = 10;              % resampling interval in ms
param.delay = 10;            % number of time delay points 
param.normType = 'norm';     % remove mean and normalize to unit variance

% Pre-processing options
param.nDim = 50;             % dimensions kept in dimension reduction (defaults to all)
param.whiten = 1;            % sets whether to whiten the data

end