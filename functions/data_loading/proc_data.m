function [data,worddata,wordboundsAN]=proc_data(y,delay,wordbounds,digitinds,digitlabels,fs,int)
% Adds delays to data and returns data divided into digits and words.
% Uses 'sliceMatrix2.m'.

% add delays
data=sliceMatrix2(y,delay); 

% get boundaries in resampled units
wordboundsAN=cell(size(wordbounds));
for i=1:length(wordbounds)
    wordboundsAN{i}=floor( (wordbounds{i}-1)/fs/int*1000+1 );
end

% slice into words and group by digits
Ndigits=length(digitlabels);
worddata=cell(1,Ndigits);
for i=1:Ndigits
    worddata{i}=cell(1,length(digitinds{i}));
    for j=1:length(worddata{i})
        thesebounds=wordboundsAN{digitinds{i}(j)};
        worddata{i}{j}=data(:,thesebounds(1):thesebounds(2));
    end
end