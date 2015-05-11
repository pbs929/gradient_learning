% Runs Gradient Descent STRF-learning experiment on Aurora-2 data.  Based 
% initially on the TrainPercep script.   
%
% Speech data mixed with noise is used for training and testing.  Receptive
% fields and objective functions are returned.  The training objectives
% are (param.costOpt):
%   -1 - no processing (control) - use with gaussian input data.
%    0 - random filters
%    1 - PCA
%    2 - ICA
%    3 - kurtosis
%    4 - sustained firing
%    5 - robustness
% 
% A control experiment can be run on gaussian data rather than speech data
% for any of the objectives (option param.gaussTest). 
%
% functions used: loadData

reload = 0;

if reload, clear all; reload=1; end
addpath(genpath('./functions'));
addpath('../software/FastICA_25')
rng shuffle; % re-seed the random number generator

%% Set parameters:

param = params_GradDesc;   % load basic options

param.gaussTest = 0;       % option to run contol experiments with random gaussian data
param.costOpt   = 1;       % *** choose the training paradigm ***
param.overwrite = 0;       % overwrite the saved file?

%% Load data

if param.gaussTest && (reload || ~exist('M','var'))    % control data (gaussian)
    [M,    MClean,    MNoise]    = loadMat_gauss(32*param.delay, 10000, param);
    [MTst, MTstClean, MTstNoise] = loadMat_gauss(32*param.delay, 10000, param);
elseif (reload || ~exist('M','var'))                   % speech data
    [M,    MClean,    MNoise] = loadMat_noise(param,'train');
    [MTst, MTstClean, MTstNoise] = loadMat_noise(param,'test');
end

%% Run training

switch param.costOpt
    
    case -1 % No processing (control)--------------------------------------
        savename = '';
        H = eye(size(M,1)); R = M; % STRF filters and response
    
    case 0  % Random filters ----------------------------------------------
        savename = 'Rand_50dim';
        [H,R] = gradDescent(M,param);
    
    case 1  % PCA ---------------------------------------------------------
        savename = 'PCA_50dim';
        [H,R] = gradDescent(M,param);
    
    case 2  % ICA ---------------------------------------------------------
        % Contrast functions: tanh, gauss, pow3
        param.contrastOpt = 'tanh';
        param.learnRate = 1;
        param.maxIter = 2000; 
        
        savename = 'ICA_50dim.mat';
        if exist(savename, 'file') && ~param.overwrite
            disp(['loading from file: ', savename])
            load(savename)
        else
            % symmetric
            [R,~,H] = fastica(M,'approach','symm','g',param.contrastOpt,'a1',1,'a2',1,'mu',param.learnRate,'maxNumIterations',param.maxIter,'lastEig',param.nDim,'epsilon',1.e-6); 
            % deflation
            %[R,~,H] = fastica(M,'approach','defl','g',param.contrastOpt,'a1',1,'a2',1,'mu',param.learnRate,'maxNumIterations',param.maxIter,'lastEig',param.nDim,'epsilon',1.e-6); 
            H = H';
        end
        
    case 3  % Kurtosis ----------------------------------------------------
        param.learnRate = 1.e1;
        param.maxIter = 1000;
        
        savename = 'Kurt_50dim_1000step.mat';
        if exist(savename, 'file') && ~param.overwrite
            disp(['loading from file: ', savename])
            load(savename)
        else
            [H,R] = gradDescent(M,param);
        end
        
    case 4  % Sustained firing --------------------------------------------
        param.deltaT = param.delay;
        param.learnRate = 1.e1;
        param.maxIter = 100; % 100 for deflation, 1000 otherwise
        
        savename = 'SF_50dim_defl.mat';
        if exist(savename, 'file') && ~param.overwrite
            disp(['loading from file: ', savename])
            load(savename)
        else
            %[H,R] = gradDescent(M,param);
            [H,R] = gradDescentDeflate(M,param);
        end
        
    case 5  % Robustness --------------------------------------------------
        param.learnRate = 1.e7;
        param.maxIter = 250;
        
        savename = 'Robust_50dim_defl.mat';
        if exist(savename, 'file') && ~param.overwrite
            disp(['loading from file: ', savename])
            load(savename)
        else
            %[H,R] = gradDescentMulti(M,MClean,MNoise,param);
            [H,R] = gradDescentDeflateMulti(M,MClean,MNoise,param);
        end
        
end

if ~isempty(savename) && (~exist(savename, 'file') || param.overwrite)
    save(savename, 'H','param'); 
end

%% get responses and show results.

% Compute testing data responses
[RTst, RTstClean, RTstNoise] = strfResponses(H, MTst, MTstClean, MTstNoise); 

% Compute cost functions
param.deltaT = param.delay;
[costSus, costKrt, costICA, costRob] = dispCosts(RTst, RTstClean, RTstNoise, param);  % display cost functions

% Display filters 
dispFilts(H, 50:-1:1, param); % select which cost by which to order. 


return

%% c++ project: export binary data (added 6-2-13)
% This exports training data as a binary file to be used with my gradient
% descent code in c++
whichFiles = [1]; % If empty, defaults to all
saveName = './trainingdata.bin';
export_binary(whichFiles, param, saveName)
