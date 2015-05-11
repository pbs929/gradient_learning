function [y,f,t,data,worddata,annot,files]=load_data_multi(sampleName,sampleDir,whichFiles,periph,normType,int,delay)
% Loads multiple spoken digit data files, processes with a peripheral 
% auditory simulation, adds delays, and slices into word-length exemplars 
% of each digit.  This is done by calling 'load_data.m' and concatenating 
% the results in a cell array.  The output format is identical to 
% 'load_data' except each variable is put into a cell array indexed by 
% speaker ID.  
% 
% inputs: 
%   sampleName - name of input samples with WILDCARDS permitted
%   sampleDir - path to input sample
%   whichFiles - specifies a subset of the directory listing - vector of indices.  
%   periph - option for peripheral simulation: spec, speccomp, speclin, gammalin, meddis, NSL
%   normType - if 'norm', normalizes each peripheral channel to unit variance; otherwise, not.  
%   int - peripheral sampling interval in milliseconds
%   delay - number of delays to include
%
% outputs: 
%   y{spkr} - peripheral model output arrayed by speaker 
%   f - frequencies in peripheral output
%   t - times in peripheral output
%   data{spkr} - peripheral output with added delays
%   worddata{spkr}{digit}{word} - peripheral output with added delays and 
%       sliced into words
%   annot{spkr} - data annotations in structure array

% get list of data files
files = rdir([sampleDir, sampleName]);
if ~isempty(whichFiles), files = files(whichFiles); end
if isempty(files), error('No input files specified'); end

% load all data into memory
y = cell(1,length(files));
data = cell(1,length(files));
worddata = cell(1,length(files));
annot = cell(1,length(files));
for fle = 1:length(files)
    disp(['Loading data file: ', files(fle).name]);
    [y{fle},f,t,data{fle},worddata{fle},annot{fle}] = ...
        load_data(files(fle).name,sampleDir,periph,normType,int,delay);
end

end
