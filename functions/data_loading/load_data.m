function [y,f,t,data,worddata,annot,sw,fs]=load_data(fileName,sampleDir,periph,normType,int,delay)
% Loads spoken digit data, processes with a peripheral auditory simulation 
% (periph_sim.m), then adds delays and slices into word-length exemplars 
% of each digit (proc_data.m).  Annotations are compiled into the structure
% array 'annot'.
% 
% inputs: 
%   fileName - name of input sample
%   sampleDir - path to input sample
%   periph - option for peripheral simulation: spec, speccomp, speclin, gammalin, meddis, NSL
%   normType - if 'norm', normalizes each peripheral channel to unit variance; otherwise, not.  
%   int - peripheral sampling interval in milliseconds
%   delay - number of delays to include
%
% outputs: 
%   y - peripheral model output
%   f - frequencies in peripheral output
%   t - times in peripheral output
%   data - peripheral output with added delays
%   worddata - peripheral output with added delays and sliced into words
%       worddata{digit}{word}
%   annot - data annotations in structure array

if isempty(findstr(filesep,fileName))
    load([sampleDir,fileName]);
else
    load(fileName)
end

% peripheral simulation
[y,f,t] = periph_sim(sw,fs,periph,int,normType);  % get peripheral response

if nargout<4, return; end

% data processing
[data,worddata,wordboundsAN] = ...     % add delays and separate into words
        proc_data(y,delay,wordbounds,dig2word,digitlabels,fs,int);
    
% output annotations
annot.dig2word = dig2word; 
annot.digitlabels = digitlabels;
annot.word2dig = word2dig;
annot.wordbounds = wordbounds;
annot.wordlabels = wordlabels;
annot.wordtimes = wordtimes;
annot.wordboundsAN = wordboundsAN;
annot.nDigits = length(dig2word);

end
