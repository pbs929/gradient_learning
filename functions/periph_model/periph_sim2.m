function [y,f,t,means,stds]=periph_sim2(sw,fs,periph,int,norm,means,stds)
% gets peripheral response as determined by the model option 'periph'
% contains one subfunction, resample2d

switch periph

    case 'spec' %conventional spectrogram, log scaling
        bandwidth=200; numChannels=32;
        [~,f,y]=spgrambw_freq(sw,fs,numChannels,'ehHdDP',bandwidth,[100 4000],60,0.001,[]); %dt=1.81/bw
        y=y';
        y=resample2d(y,round(int/1000/0.001));
        t=(0:size(y,2)-1)*int/1000;   
    case 'speccomp50' %conventional spectrogram, square root compression
        bandwidth=50; numChannels=32;
        [~,f,y]=spgrambw_freq(sw,fs,numChannels,'ehHP',bandwidth,[100 4000],[],0.001,[]); %dt=1.81/bw
        y=y'.^0.25;
        y=resample2d(y,round(int/1000/0.001));
        t=(0:size(y,2)-1)*int/1000;
    case 'speccomp100' %conventional spectrogram, square root compression
        bandwidth=100; numChannels=32;
        [~,f,y]=spgrambw_freq(sw,fs,numChannels,'ehHP',bandwidth,[100 4000],[],0.001,[]); %dt=1.81/bw
        y=y'.^0.25;
        y=resample2d(y,round(int/1000/0.001));
        t=(0:size(y,2)-1)*int/1000;
    case 'speccomp200' %conventional spectrogram, square root compression
        bandwidth=200; numChannels=32;
        [~,f,y]=spgrambw_freq(sw,fs,numChannels,'ehHP',bandwidth,[100 4000],[],0.001,[]); %dt=1.81/bw
        y=y'.^0.25;
        y=resample2d(y,round(int/1000/0.001));
        t=(0:size(y,2)-1)*int/1000;
    case 'speccomp400' %conventional spectrogram, square root compression
        bandwidth=400; numChannels=32;
        [~,f,y]=spgrambw_freq(sw,fs,numChannels,'ehHP',bandwidth,[100 4000],[],0.001,[]); %dt=1.81/bw
        y=y'.^0.25;
        y=resample2d(y,round(int/1000/0.001));
        t=(0:size(y,2)-1)*int/1000;    
    case 'speclin' %conventional spectrogram, linear scaling  *** not updated
        bandwidth=200; numChannels=32;
        [t,f,y]=spgrambw_freq(sw,fs,numChannels,'ehHP',bandwidth,[100 4000],[],int/1000,[]); %dt=1.81/bw
        y=y'.^0.5;
    case 'specstdcomp' %conventional spectrogram, square root compression
        bandwidth=200; numChannels=32;
        [~,f,y]=spgrambw_freq(sw,fs,numChannels,'hHP',bandwidth,[100 4000],[],0.001,[]); %dt=1.81/bw
        y=y'.^0.25;
        y=resample2d(y,round(int/1000/0.001));
        t=(0:size(y,2)-1)*int/1000;    
    case 'specstd' %conventional spectrogram, linear frequency axis, log scaling.  
        bandwidth=200; numChannels=32;
        [~,f,y]=spgrambw_freq(sw,fs,numChannels,'hHdDP',bandwidth,[100 4000],60,0.001,[]); %dt=1.81/bw
        y=y';
        y=resample2d(y,round(int/1000/0.001));
        t=(0:size(y,2)-1)*int/1000;    
        
    case 'gamma' %gammatone filter, square root compression
        normfactor=1; numChannels=32;
        lowFreq=100; highFreq=4000;
        [y, f] =combinedANF(sw*normfactor, fs, numChannels, highFreq, lowFreq);
        y=max(y,0); y=y.^0.5;
        y=resample2d(y,round(fs*int/1000));
        t=(0:size(y,2)-1)*int/1000;    
    case 'gammalin' %gammatone, no compression
        normfactor=1; numChannels=32;
        lowFreq=100; highFreq=4000;
        [y, f] =combinedANF(sw*normfactor, fs, numChannels, highFreq, lowFreq);
        y=max(y,0); 
        y=resample2d(y,round(fs*int/1000));
        t=(0:size(y,2)-1)*int/1000;
        
    case 'meddis' %meddis model
        normfactor=1; numChannels=32;
        lowFreq=100; highFreq=4000;
        [~, f, y] =combinedANF(sw*normfactor, fs, numChannels, highFreq, lowFreq);
        y=resample2d(y,round(fs*int/1000));
        t=(0:size(y,2)-1)*int/1000;
        
    case 'NSL' %auditory spectrogram
        sw=resample(sw,2,5); fs=fs*2/5; % decrease sampling rate of input
        [y, f, int] =combinedNSL(sw, fs, int);
        y=log(y); y=max(y,max(y(:))-7); 
        y=y(:,2:2:end); f=f(2:2:end); y=y'; %y=resample2d(y',2,1); f=decimate(f,2);
        
    otherwise 
        error('invalid choice of peripheral model')
end

if nargin==5 && strcmp(norm,'norm')
    disp('periph_sim: peripheral simulation with unfixed means/variances');
    means=mean(y,2);
    stds=std(y,0,2);
    y=(y-means*ones(1,size(y,2)))./(stds*ones(1,size(y,2)));
end

if nargin>5 && strcmp(norm,'norm')
    disp('periph_sim: peripheral simulation with fixed means/variances');
    y=(y-means*ones(1,size(y,2)))./(stds*ones(1,size(y,2)));
end

if nargin<5
    disp('periph_sim: peripheral simulation without normalization')
    means=[]; stds=[];
end

end

function M2=resample2d(M,xint)
% resamples 2d matrix along rows

M2=(resample(M',1,xint))';

end
