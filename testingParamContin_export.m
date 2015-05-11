function param = testingParamContin_export()

% Directory for saving temporary data:
param.tempDir = './tempsave_cont/';

% Feature detector parameters - Set which of the feature detectors stored
% from training to use in the recognition.  The 'percepTag' determines
% which training run to used.  The 'percepDir' is the location where the
% feature detector parameters are stored.  The number of male and female 
% speakers determines which of the files to load.  
param.percepTag = 'gamma_indiv_8ms_8del_1pos_1pt_5neg_C_multi';%'gamma_norm_8ms_8del_1pos_1pt_5neg_C_multi_A'; 
param.percepDir = './trained_models/';
param.nMalePercep = 50;
param.nFemalePercep = 50;

% Template parameters - Give the path to the AURORA training data and set
% how many exemplars of each digit from each speaker to use as templates.
param.trainDir = '../data/samples_aurora/aurora_samps_train/';
param.nMale = -1;         % number of male SAMPLES to use, -1 defaults to all
param.nFemale = -1;       % number of female SAMPLES to use, -1 defaults to all

% Test data parameters - Give the path to the AURORA testing data and set
% how many exemplars of each digit from each speaker to use as testing 
% utterances.  
param.testDir = '../data/samples_aurora/aurora_samps_test/';
param.nMaleTest = -1;     % number of male SAMPLES to use, -1 defaults to all
param.nFemaleTest = -1;   % number of female SAMPLES to use -1 defaults to all

% Noise parameters - Set which noise types and levels to use for testing
param.noiseType = [1:4];    % 4 types total
param.noiseLevs = 20:-5:-5;
param.runClean = 1;       % sets whether to run clean data