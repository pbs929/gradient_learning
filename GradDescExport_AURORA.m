% Takes STRF training results from Grad Descent program and exports the
% neural responses to speech data for use in the Aurora recognition
% experiments. 

clear all
addpath(genpath('./functions'))
addpath(genpath('./Testing'))
rng shuffle; % re-seed the random number generator

%% Export sequence code as HTK features

% We save the data in HTK format for use in the paper.  A different
% directory should be used for any given parameter set.

param = testingParamContin_export();

% STRF file:
strfs = 'Robust_50dim_defl_half';

% export directory:
dirname = ['../data/samples_aurora/aurora_samps_htk_',strfs,'/']; % 50 speakers

export_htk_contin(strfs, param, dirname);
