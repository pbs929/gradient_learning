function [t,f,b]=spgrambw_freq(s,fs,nfrq,varargin)
% altered to have nfrq as an input argument PBS 11-4-12
%
%SPGRAMBW Draw spectrogram [T,F,B]=(s,fs,mode,bw,fmax,db,tinc,ann)
%
%  Usage: spgrambw(s,fs,'pJcw')  % Plot spectrogram with my favourite set of options
%         
%         For examples of the many options available see: 
%         http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/tutorial/spgrambw/spgram_tut.pdf
%
%  Inputs:  S         speech signal, or single-sided power spectrum array, S(NT,NF), in power per Hz
%           FS        sample fequency (Hz) or [FS T1] where T1 is the time of the first sample
%                     or, if s is a matrix, [FS T1 FINC F1] where FS is the frame rate, T1 is
%                     the time of the first sample, FINC is the frequency increment and F1 the
%                     frequency of the first column.
%           MODE      optional character string specifying options (see list below)
%           BW        bandwidth resolution in Hz (DFT window length = 1.81/BW)[default: 200]
%           FMAX      frequency range [Fmin Fstep Fmax]. If Fstep is omitted
%                     it is taken to be (Fmax-Fmin)/257, if Fmin is also omitted it is taken
%                     to be 0 (or 20Hz for mode l), if all three are omitted Fmax is taken to be FS/2.
%                     If modes m, b, e or l are specified then the units are in mel, bark or erb or
%                     log10(Hz); this can be over-ridden by the 'h' option.
%           DB        either dB-range or [dB-min dB-max] [default: 40]
%           TINC      output frame increment in seconds [0 or missing uses default=0.45/BW]
%                     or [TFIRST TLAST] or [TFIRST TINC TLAST] where TFIRST/TLAST are the times
%                     of first/last frames
%           ANN       annotation cell array: each row contains either
%                     {time 'text-string' 'font'} or {[t_start t_end] 'text-string' 'font'} where
%                     the time value is in seconds with s(n) at time offset+n/fs. The font column can
%                     omitted in which case the system font will be used. MATLAB cannot cope with
%                     unicode so I recommend the SILDoulosIPA (serifed) or SILSophiaIPA (sans) fonts
%                     for phonetic symbols; these are now a little hard to find.
%
% Outputs:  T(NT)        time axis values (in seconds). Input sample s(n) is at time offset+n/fs.
%           F(NF)        frequency axis values in Hz or, unless mode=H, other selected frequency units
%                        according to mode: m=mel, l=log10(Hz), b=bark,e=erb-rate
%           B(NT,NF)     spectrogram values in power (or clipped dB values if 'd' option given)
%
% MODE:  'p' = output power per decade rather than power per Hz [preemphasis]
%        'P' = output power per mel/bark/erb according to y axis scaling
%        'd' = output B array is in dB rather than power
%        'D' = clip the output B array to the limits specified by the "db" input
%
%        'm' = mel scale
%        'b' = bark scale
%        'e' = erb scale
%        'l' = log10 Hz frequency scale
%        'f' = label frequency axis in Hz rather than mel/bark/... 
%
%        'h' = units of the FMAX input are in Hz instead of mel/bark
%              [in this case, the Fstep parameter is used only to determine
%               the number of filters]
%        'H' = express the F output in Hz instead of mel/bark/...
%
%        'g' = draw a graph even if output arguments are present
%        'j' = jet colourmap
%        'J' = "thermal" colourmap that is linear in grayscale. Based on Oliver Woodford's
%                 real2rgb at http://www.mathworks.com/matlabcentral/fileexchange/23342
%        'i' = inverted colourmap (white background)
%        'c' = include a colourbar as an intensity scale
%        'w' = draw the speech waveform above the spectrogram
%        'a' = centre-align annotations rather than left-aligning them
%        't' = add time markers with annotations
%
% The BW input gives the 6dB bandwidth of the Hamming window used in the analysis.
% Equal amplitude frequency components are guaranteed to give separate peaks if they
% are this far apart. This value also determines the time resolution: the window length is
% 1.81/BW and the low-pass filter applied to amplitude modulations has a 6-dB bandwidth of
% BW/2 Hz.
%
% The units are power per Hz unless the u
% option is given in which case power per displayed unit is used
% or power per decade for the l option.

%%%% BUGS %%%%%%
% * allow ANN rows to be a mixture of intervals and instants
% * allow multiple ANN rows
% * Do not use triangular interpolation if the output frequencies are the same as an FFT
% * Place as many subticks as will fit beyond the last tick with the 'f' option
% * Use a special subtick pattern between ticks that are powers of 10 using the 'f' option
% * Future options:
%       ['q' = constant q transform]
%       ['k' = add a piano keyboard to the frequency scale]
%       ['z' = use a bipolar colourmap for a matrix input with negative values]

%      Copyright (C) Mike Brookes 1997-2011
%      Version: $Id: spgrambw.m,v 1.22 2011/06/07 21:00:55 dmb Exp $
%
%   VOICEBOX is a MATLAB toolbox for speech processing.
%   Home page: http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You can obtain a copy of the GNU General Public License from
%   http://www.gnu.org/copyleft/gpl.html or by writing to
%   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
persistent tcmap
if isempty(tcmap)
    % modified thermal with better grayscale linearity
    tcmap=[ 0 0 0; 7 0 17; 14 0 33; 21 0 50; 29 0 67; 36 0 84; 43 0 100; 50 0 117;
        57 0 134; 64 0 150; 72 0 167; 80 3 164; 89 7 156; 97 11 149; 106 15 142; 114 19 134;
        123 23 127; 131 27 119; 140 31 112; 149 35 105; 157 39 97; 166 43 90; 174 47 82;
        183 51 75; 192 55 68; 200 59 60; 209 63 53; 217 67 45; 226 71 38; 234 75 31;
        243 79 23; 252 83 16; 255 88 12; 255 95 12; 255 102 11; 255 109 11; 255 116 10;
        255 123 10; 255 130 9; 255 137 9; 255 144 8; 255 151 8; 255 158 7; 255 165 7;
        255 172 6; 255 179 6; 255 186 5; 255 193 4; 255 200 4; 255 207 3; 255 214 3; 255 221 2;
        255 228 2; 255 235 1; 255 242 1; 255 249 0; 255 252 22; 255 252 55; 255 253 88;
        255 253 122; 255 254 155; 255 254 188; 255 255 222; 255 255 255]/255;
end
if nargin<2
    error('Usage: SPGRAMBW(s,fs,mode,bw,fmax,db,tinc)');
end
%SPGRAMBW Draw grey-scale spectrogram [T,F,B]=(s,fs,mode,bw,fmax,db,tinc)
%
% first decode the input arguments
%
if size(s,1)==1
    s=s(:);   % force to be a column vector (unless it is a matrix)
end
[ns1,ns2]=size(s);
ap=zeros(1,6);
j=2;
if numel(fs)<2
    fs(2)=1/fs(1);  % first sample or frame is at time 1/fs
end
for i=1:length(varargin)
    if ischar(varargin{i})
        ap(1)=i;
    else
        ap(j)=i;
        j=j+1;
    end
end
if ap(1) && ~isempty(varargin{ap(1)})
    mode=varargin{ap(1)};
else
    mode='';  % default mode
end
if ap(2) && ~isempty(varargin{ap(2)})
    bw=varargin{ap(2)};
else
    bw=200;
end
if ap(3) && ~isempty(varargin{ap(3)})
    fmax=varargin{ap(3)};
else
    fmax=[];
end
if ap(4) && ~isempty(varargin{ap(4)})
    db=varargin{ap(4)};
else
    db=40;
end
if ap(5) && ~isempty(varargin{ap(5)})
    tinc=varargin{ap(5)};
else
    tinc=0;
end
switch numel(tinc)
    case 1
        tinc=[tinc -Inf Inf];
    case 2
        tinc=[0 tinc];
    otherwise
        tinc=tinc([2 1 3]);
end
if tinc(1)<=0
    tinc(1)=1.81/(4*bw); % default frame increment
end
if ap(6)
    ann=varargin{ap(6)};
else
    ann=[];
end

% now sort out the mode flags

mdsw='  ';           % [yscale preemph]
for i=1:length(mode)
    switch mode(i)
        case {'l','m','b','e'}
            mdsw(1)=mode(i);
        case {'p','P'}
            mdsw(2)=mode(i);
    end
end
if mdsw(2)=='P'
    mdsw(2)=mdsw(1);        % preemphasis is scaling dependent
end
%
% sort out the frequency axis
%
flmin=30;                   % min frequency for 'l' option
if ns2==1
    fnyq=fs(1)/2;           % default upper frequency limit is fs/2
else                        % input is a power spectrum
    if numel(fs)<3
        fs(3)=fs(1)*0.25;   % default increment is 0.25 times frame increment
    end
    if numel(fs)<4
        fs(4)=0;            % first freq bin is DC by default
    end
    fnyq=fs(4)+(ns2-1)*fs(3);  % default upper frequency limit is highest supplied frequency
end

if ~numel(fmax)             % no explicit frequency range
    switch mdsw(1)
        case 'l'
            fx=linspace(log10(flmin),log10(fnyq),nfrq);   % 20  Hz to Nyquist
        case 'm'
            fx=linspace(0,frq2mel(fnyq),nfrq);   % DC to Nyquist
        case 'b'
            fx=linspace(0,frq2bark(fnyq),nfrq);   % DC to Nyquist
        case 'e'
            fx=linspace(0,frq2erb(fnyq),nfrq);   % DC to Nyquist
        otherwise   % linear Hz scale
            fx=(0:nfrq-1)*fnyq/(nfrq-1);
    end
else
    if any(mode=='h')
        switch mdsw(1)
            case 'l'
                fmaxu=log10(fmax);   % 20  Hz to Nyquist
            case 'm'
                fmaxu=frq2mel(fmax);   % DC to Nyquist
            case 'b'
                fmaxu=frq2bark(fmax);   % DC to Nyquist
            case 'e'
                fmaxu=frq2erb(fmax);   % DC to Nyquist
            otherwise
                fmaxu=fmax;  % linear Hz scale
        end
    else
        fmaxu=fmax;                 % already in the correct units
    end
    if numel(fmax)<2   % only max value specified
        if mdsw(1)=='l'
            fx=linspace(log10(flmin),fmaxu,nfrq);   % 20  Hz to fmax
        else
            fx=linspace(0,fmaxu,nfrq);   % DC to fmax
        end
    elseif numel(fmax)<3 % min and max values specified
        fx=linspace(fmaxu(1),fmaxu(2),nfrq);   % fmin to fmax
    else
        fmaxu(2)=fmax(2)*(fmaxu(3)-fmaxu(1))/(fmax(3)-fmax(1)); % scale the step size appropriately
        fx=fmaxu(1):fmaxu(2):fmaxu(3);   % fmin to fmax in steps of finc
        nfrq=length(fx);
    end
end
switch mdsw(1)          % convert the frequency range to Hz
    case 'l'
        f=10.^fx;
        frlab='log_{10}Hz';
        frlabf='log';
        frq2y=@log10;
        y2frq=@(x) 10.^x;
    case 'm'
        f=mel2frq(fx);
        frlab='Mel';
        frlabf='Mel';
        frq2y=@frq2mel;
        y2frq=@mel2frq;
    case 'b'
        f=bark2frq(fx);
        frlab='Bark';
        frlabf='Bark';
                frq2y=@frq2bark;
        y2frq=@bark2frq;
    case 'e'
        f=erb2frq(fx);
        frlab='Erb-rate';
        frlabf='Erb';
        frq2y=@frq2erb;
        y2frq=@erb2frq;
    otherwise
        f=fx;
        frlab='Hz';
                frq2y=@(x) x;
        y2frq=@(x) x;
end
if ~any(mode=='H')
    f=fx;               % give output frequencies in native units instead of Hz unless 'H' is specified
end
%
% now calculate the spectrogram
%
if ns2==1   % input is a speech signal vector
    winlen = fix(1.81*fs(1)/bw);   % window length
    win=0.54+0.46*cos((1-winlen:2:winlen)*pi/winlen);  % Hamming window
    ninc=max(round(tinc(1)*fs(1)),1);                 % window increment in samples
    %  we need to take account of minimum freq increment + make it exact if possible
    fftlen=pow2(nextpow2(4*winlen));        % enough oversampling to get good interpolation
    win=win/sqrt(sum(win.^2));              % ensure window squared sums to unity
    ix1=max(round((tinc(2)-fs(2))*fs(1)-(winlen-3)/2),1); % first sample required
    ix2=min(ceil((tinc(3)-fs(2))*fs(1)+(winlen+1)/2),ns1); % last sample required
    [sf,t]=enframe(s(ix1:ix2),win,ninc);
    t=fs(2)+(t+ix1-2)/fs(1);                         % time axis
    b=rfft(sf,fftlen,2);
    b=b.*conj(b)*2/fs(1);          % Power per Hz
    b(:,1)=b(:,1)*0.5;   % correct for no negative zero frequency to double the power
    b(:,end)=b(:,end)*0.5;   % correct for no negative nyquist frequency to double the power
    fb=(0:fftlen/2)*fs(1)/fftlen; % fft bin frequencies
    fftfs=fs(1);
else

    b=s;
    t=fs(2)+(0:ns1-1)/fs(1);  % frame times
    fb=fs(4)+(0:ns2-1)*fs(3);
    fftlen=[ns2 fs(3) fs(4)]; % for filtbankm: ns2=# input freq bins, freq increment (Hz), first bin freq (Hz)
    fftfs=0;
    %     fftlen=2*(ns2-1);  % assume an even length fft
    %     fftfs=fftlen*fs(3);
end
nfr=numel(t);                   % number of frames
dblab='Power/Hz';
switch mdsw(2)
    case {'p','l'}
        b=b.*repmat(fb*log(10),nfr,1);       % convert to power per decade
        dblab='Power/Decade';
    case 'm'
        b=b.*repmat((1+fb/700)*log(1+1000/700)/1000,nfr,1);       % convert to power per mel
        dblab='Power/Mel';
    case 'b'
        b=b.*repmat((1960+fb).^2/52547.6,nfr,1);       % convert to power per bark
        dblab='Power/Bark';
    case 'e'
        b=b.*repmat(6.23*fb.^2 + 93.39*fb + 28.52,nfr,1);       % convert to power per erb
        dblab='Power/Erb-rate';
end
%
% Now map onto the desired frequency scale
%
b=b*filtbankm(nfrq,fftlen,fftfs,fx(1),fx(end),['cush' mdsw(1)])';

if ~nargout || any(mode=='g') ||  any(mode=='d')
    if numel(db)<2          % find clipping limits
        plim=max(b(:))*[0.1^(0.1*db) 1];
    else
        plim=10.^(0.1*db(1:2));
    end
    if plim(2)<=0
        plim(2)=1;
    end
    if plim(1)<=0 || plim(1)==plim(2)
        plim(1)=0.1*plim(2);
    end
    if ~nargout || any(mode=='g')
        bd=10*log10(b);  % save an unclipped log version for plotting
    end
    if any(mode=='D')
        b=min(max(b,plim(1)),plim(2)); % clip the output
    end
    if any(mode=='d')
        b=10*log10(b);    % output the dB version
    end
end
% now plot things
if ~nargout || any(mode=='g')
    cla;  % clear current axis
    imagesc(t,fx,bd');
    axis('xy');
    set(gca,'tickdir','out','clim',10*log10(plim));
    if any(mode=='j')
        colormap('jet');
        map=colormap;
    elseif any(mode=='J')
        map=tcmap;
    else
        map = repmat((0:63)'/63,1,3);
    end
    if any(mode=='i')               % 'i' option = invert the colourmap
        map=map(64:-1:1,:);
    end
    colormap(map);
    if any(mode=='c')                % 'c' option = show a colourbar
        colorbar;
        cblabel([dblab ' (dB)']);
    end
    %
    % Now check if annotations or a waveform are required
    %
    dotaw=[((any(mode=='t') && size(ann,2)>1) || size(ann,2)==1) size(ann,2)>1 (any(mode=='w') && ns2==1)];
        ylim=get(gca,'ylim');
    if  any(dotaw)
        yrange = ylim(2)-ylim(1);
        zlim=ylim;
        toptaw=cumsum([0 dotaw.*[0.05 0.05 0.1]]*yrange)+ylim(2);
        zlim(2)=toptaw(4);
        set(gca,'ylim',zlim,'color',map(1,:));
        if dotaw(3)        % Plot the waveform
            smax=max(s(:));
            smin=min(s(:));
            srange=smax-smin;
            hold on
            plot(fs(2)+(0:length(s)-1)/fs(1),(s-smin)/srange*0.9*(toptaw(4)-toptaw(3))+toptaw(3),'color',map(48,:))
            hold off
        end
        if dotaw(1) || dotaw(2)
            tmk=cell2mat(ann(:,1));
            tmksel=tmk(:,1)<=t(end) & tmk(:,end)>=t(1);
            yix=1+[tmk(tmksel,1)<t(1) ones(sum(tmksel),2) tmk(tmksel,end)>t(end)]';
            tmk(:,1)=max(tmk(:,1),t(1));  % clip to axis limits
            tmk(:,end)=min(tmk(:,end),t(end));
        end
        if dotaw(1) && any(tmksel)  % draw time markers
            ymk=toptaw(1:2)*[0.8 0.4;0.2 0.6];
            switch size(tmk,2)
                case 0
                case 1      % isolated marks
                    hold on
                    plot([tmk(tmksel) tmk(tmksel)]',repmat(ymk',1,sum(tmksel)),'color',map(48,:));
                    hold off
                otherwise % draw durations

                    hold on
                    plot(tmk(tmksel,[1 1 2 2])',ymk(yix),'color',map(48,:));
                    hold off
            end
        end
        if dotaw(2) && any(tmksel) % print annotations
            if any(mode=='a')
                horal='center';
                tmk=(tmk(:,1)+tmk(:,end))*0.5;
            else
                horal='left';
                tmk=tmk(:,1);
            end
            if size(ann,2)>2
                font='Arial';
                for i=1:size(ann,1)
                    if tmksel(i)
                        if ~isempty(ann{i,3})
                            font = ann{i,3};
                        end
                        text(tmk(i),toptaw(2),ann{i,2},'color',map(48,:),'fontname',font,'VerticalAlignment','baseline','HorizontalAlignment',horal);
                    end
                end
            else
                for i=1:size(ann,1)
                    if tmksel(i)
                        text(tmk(i),toptaw(2),ann{i,2},'color',map(48,:),'VerticalAlignment','baseline','HorizontalAlignment',horal);
                    end
                end
            end
        end
    end
    xlabel(['Time (' xticksi 's)']);
    if any(mode=='f') && ~strcmp(frlab,'Hz')
        ylabel([frlabf '-scaled frequency (Hz)']);
        ytickhz(frq2y,y2frq);
    else
    ylabel(['Frequency (' yticksi frlab ')']);
    end
    ytick=get(gca,'YTick');
    ytickl=get(gca,'YTickLabel');
    msk=ytick<=ylim(2);
    if any(~msk)
        set(gca,'YTick',ytick(msk),'YTickLabel',ytickl(msk));
    end
end

function ytickhz(frq2y,y2frq)
% label non linear y frequency axis
%
% Bugs/Suggestions:
% * Add a penalty for large numbers (e.g. 94 is less "round" than 11)
% * possibly add subticks at 1:2:5 if boundaries are 1 and 10
% * could treat subtick allocation specially if bounding lables are both powers of 10
%   and work in log spacing rather than spacing directly

% algorithm constants

seps=[0.4 1 3 6]; % spacings: (a) min subtick, (b) min tick, (c) min good tick, (d) max good tick
ww=[0.5 0.6 0.8 0.1 0.3 0.3 0.2];  % weight for (a) last digit=5, (b) power of 10, (c) power of 1000, (d) equal spacing, (e) 1:2:5 labels (f) <seps(3) (g) >seps(4)
nbest=10; % number of possibilities to track

prefix={'y','z','a','f','p','n','�','m','','k','M','G','T','P','E','Z','Y'};

ah=gca;
getgca=get(ah);  % Get original axis properties
set(ah,'Units','points','FontUnits','points');
getgcac=get(ah);  % Get axis properties in points units
set(ah,'Units',getgca.Units,'FontUnits',getgca.FontUnits); % return to original values
ylim=getgca.YLim;
yrange=ylim*[-1;1];
chsz= yrange*getgcac.FontSize/getgcac.Position(4); % char height in Y-units
% divide the y-axis up into bins containing at most one label each
maxl=ceil(2*yrange/chsz);  % max number of labels

% candidate array [cand(:,[1 2])/1000 cand(:,5) cand(:,6)/1000 cand(:,[7 8])]
% 1,2=y limits, 3,4=log limits, 5=Hz, 6=cost, 7=mantissa, 8=exponent, 9=sig digits, 10=y-position
cand=zeros(maxl+2,10);
yinc=(yrange+chsz*0.0002)/maxl;  % bin spacing (allowing for a tiny bit to ensure the ends are included)
cand(2:end-1,2)=ylim(1)+yinc*(1:maxl)'-chsz*0.0001;
cand(3:end-1,1)=cand(2:end-2,2);
cand(2,1)=cand(2,2)-yinc;
cand(2:end-1,1:2)=y2frq(max(cand(2:end-1,1:2),0));

% find the "roundest" number in each interval
% first deal with intervals containing zero
cand([1 maxl+2],6)=-1;
cand(2,9)=(cand(2,1)<=0);  % mask out interval contaiing zero
cand(2,6)=-cand(2,9);
msk=cand(:,6)==0;  % find rows without a cost yet
cand(msk,3:4)=log10(cand(msk,1:2));
% find powers of 1000
loglim=ceil(cand(:,3:4)/3);
msk=loglim(:,2)>loglim(:,1);
if any(msk)
    xp=loglim(msk,1);
    wuns=ones(length(xp),1);
    cand(msk,5:9)=[1000.^xp wuns-ww(3) wuns 3*xp wuns];
end
% find powers of 10
loglim=ceil(cand(:,3:4));
msk=~msk & (loglim(:,2)>loglim(:,1));
if any(msk)
    xp=loglim(msk,1);
    wuns=ones(length(xp),1);
    cand(msk,5:9)=[10.^xp wuns-ww(2) wuns xp wuns];
end
% find value with fewest digits
msk=cand(:,6)==0;  % find rows without a cost yet
maxsig=1-floor(log10(10^min(cand(msk,3:4)*[-1;1])-1)); % maximum number of significant figures to consider
pten=10.^(0:maxsig-1);   % row vector of powers of ten
noten=10.^(-floor(cand(msk,3))); % exponent of floating point representation of lower bound
sigdig=sum((ceil(cand(msk,2).*noten*pten)-ceil(cand(msk,1).*noten*pten))==0,2); % number of digits common to the interval bounds
lowman=ceil(cand(msk,1).*noten.*10.^sigdig);
midman=10*floor(lowman/10)+5;
highman=ceil(cand(msk,2).*noten.*10.^sigdig);
mskman=midman>=lowman & midman<highman;   % check if we can include a manitssa ending in 5
lowman(mskman)=midman(mskman);
cand(msk,6:9)=[sigdig+1 lowman floor(cand(msk,3))-sigdig sigdig+1];
cand(msk,5)=cand(msk,7).*10.^cand(msk,8);
cand(msk,6)=cand(msk,6)-(mod(cand(msk,7),10)==5)*ww(1);
cand(2:end-1,10)=frq2y(cand(2:end-1,5));
cand([1 maxl+2],10)=ylim + seps(4)*chsz*[-1 1]; % put imaginary labels at the optimum spacing beyond the axes
% [cand(:,[1 2 5])/1000 cand(:,[6 7 8 9])]

% Now do n-best DP to find the best sequence

ratint=[8/5 25/10 0 0 4/3];
costs=repmat(Inf,nbest,maxl+2); % cumulative path costs
costs(1,1)=0; % starting node only has one option
prev=ones(nbest,maxl+2); % previous label in path
labcnt=zeros(nbest,maxl+2); % number of labels in path
for i=2:maxl+2
    ntry=nbest*(i-1); % number of previous options
    prevc=reshape(repmat(1:i-1,nbest,1),ntry,1); % previous candidate
    prevprev=1+floor((prev(1:ntry)'-1)/nbest); % previous previous candidate
    msk=prevprev>1+(maxl+2)*(i==maxl+2); % mask for label triplets
    labcnti=labcnt(1:ntry)+1;
    disti=(cand(i,10)-cand(prevc,10))/chsz; % distance to previous label in characters
    costa=max(seps(3)-disti,0)*ww(6)+max(disti-seps(4),0)*ww(7);
    incri=(cand(i,5)-cand(prevc,5)); % label increment
    incrj=(cand(i,5)-cand(prevprev,5)); % double label increment
    if any(msk)
        costa(msk)=costa(msk)- ww(4)*(abs(incrj(msk)-2*incri(msk))<0.01*incri(msk));
        if cand(i,7)==1 || cand(i,7)==2 || cand(i,7)==5 % look for labels 1:2:5
            costa(msk)=costa(msk)- ww(5)*(abs(incrj(msk)-ratint(cand(i,7))*incri(msk))<0.01*incri(msk));
        end
    end
    costa(disti<seps(2))=Inf;
    costi=(costs(1:ntry).*max(labcnt(1:ntry),1)+costa'+cand(i,6))./labcnti;
    [sc,isc]=sort(costi);
    isc=isc(1:nbest);
    costs(:,i)=sc(1:nbest)';
    prev(:,i)=isc';
    labcnt(:,i)=labcnti(isc)';
end

% now traceback the best sequence

% fprintf('Traceback\n\n');
ichoose=0;
labchoose=[];
for i=1:nbest
    if labcnt(i,maxl+2)>1 && costs(i,maxl+2)<Inf
        lablist=zeros(labcnt(i,maxl+2)-1,1);
        k=prev(i,maxl+2);
        for j=labcnt(i,maxl+2)-1:-1:1
            lablist(j)=1+floor((k-1)/nbest);
            k=prev(k);
        end
%         fprintf('Cost=%8.2f :',costs(i,maxl+2));
%         fprintf(' %g',cand(lablist,5))
%         fprintf('\n');
        if ~ichoose || labcnt(ichoose,maxl+2)==1
            ichoose=i;
            labchoose=lablist;
        end
    end
end

% now create the labels

ntick=length(labchoose);
% sort out the subticks
subpos=[];
if ntick>=2
    for i=1:ntick-1
        clj=cand(labchoose(i:i+1),:);
        sprec=min(clj(1,8)+100*(clj(1,7)==0),clj(2,8)); % subtick precision
        spos=(clj(1,7)*10^(clj(1,8)-sprec):clj(2,7)*10^(clj(2,8)-sprec))*10^sprec;
        nsub=length(spos);
        if nsub==2
            spos=spos*[1 0.5 0;0 0.5 1];
            nsub=3;
        end
        if nsub>=3
            yspos=frq2y(spos);
            for kk=1:3 % try various subdivisions: every 1, 2 or 5
                k=kk+2*(kk==3);  % 1, 2 and 5
                if 2*k<=nsub-1 && ~mod(nsub-1,k)  % must divide exactly into nsub
                    if all((yspos(1+k:k:nsub)-yspos(1:k:nsub-k))>=(seps(1)*chsz)) % check they all fit in
                        subpos=[subpos yspos(1+k:k:nsub-k)];
                        if i==1
                            spos=(ceil(cand(2,1)/10^sprec):clj(1,7)*10^(clj(1,8)-sprec))*10^sprec;
                            nsub=length(spos);
                            yspos=frq2y(spos);
                            if nsub>=k+1 && all((yspos(nsub:-k:1+k)-yspos(nsub-k:-k:1))>=(seps(1)*chsz))
                                subpos=[subpos yspos(nsub-k:-k:1)];
                            end
                        elseif i==ntick-1
                            spos=(clj(2,7)*10^(clj(2,8)-sprec):floor(cand(end-1,2)/10^sprec))*10^sprec;
                            nsub=length(spos);
                            yspos=frq2y(spos);
                            if nsub>=k+1 && all((yspos(1+k:k:nsub)-yspos(1:k:nsub-k))>=(seps(1)*chsz))
                                subpos=[subpos yspos(1+k:k:nsub)];
                            end
                        end
                        break;
                    end
                end
            end
        end
    end
end
nsub=length(subpos);
tickpos=[cand(labchoose,10); subpos'];
ticklab=cell(ntick+nsub,1);
sipref=min(max(floor((sum(cand(labchoose,8:9),2)-1)/3),-8),8);
nzadd=cand(labchoose,8)-3*sipref;  % trailing zeros to add
digzer=cand(labchoose,7).*10.^max(nzadd,0); % label digits including trailing zeros
ndleft=cand(labchoose,9)+nzadd; % digits to the left of the decimal point
for i=1:ntick
    tickint=num2str(digzer(i));
    if nzadd(i)<0
        tickint=[tickint(1:ndleft(i)) '.' tickint(1+ndleft(i):end)];
    end
    ticklab{i} = sprintf('%s%s',tickint,prefix{sipref(i)+9});
end
for i=ntick+1:ntick+nsub
    ticklab{i}='';
end
[tickpos,ix]=sort(tickpos);
ticklab=ticklab(ix);

set(ah,'YTick',tickpos','YTickLabel',ticklab);

function [mel,mr] = frq2mel(frq)
%FRQ2ERB  Convert Hertz to Mel frequency scale MEL=(FRQ)
%	[mel,mr] = frq2mel(frq) converts a vector of frequencies (in Hz)
%	to the corresponding values on the Mel scale which corresponds
%	to the perceived pitch of a tone.
%   mr gives the corresponding gradients in Hz/mel.

%	The relationship between mel and frq is given by:
%
%	m = ln(1 + f/700) * 1000 / ln(1+1000/700)
%
%  	This means that m(1000) = 1000
%
%	References:
%
%	  [1] S. S. Stevens & J. Volkman "The relation of pitch to
%		frequency", American J of Psychology, V 53, p329 1940
%	  [2] C. G. M. Fant, "Acoustic description & classification
%		of phonetic units", Ericsson Tchnics, No 1 1959
%		(reprinted in "Speech Sounds & Features", MIT Press 1973)
%	  [3] S. B. Davis & P. Mermelstein, "Comparison of parametric
%		representations for monosyllabic word recognition in
%		continuously spoken sentences", IEEE ASSP, V 28,
%		pp 357-366 Aug 1980
%	  [4] J. R. Deller Jr, J. G. Proakis, J. H. L. Hansen,
%		"Discrete-Time Processing of Speech Signals", p380,
%		Macmillan 1993
%	  [5] HTK Reference Manual p73
%	



%      Copyright (C) Mike Brookes 1998
%      Version: $Id: frq2mel.m,v 1.7 2010/08/01 08:41:57 dmb Exp $
%
%   VOICEBOX is a MATLAB toolbox for speech processing.
%   Home page: http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You can obtain a copy of the GNU General Public License from
%   http://www.gnu.org/copyleft/gpl.html or by writing to
%   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
persistent k
if isempty(k)
    k=1127.01048;
end
af=abs(frq);
mel = sign(frq).*log(1+af/700)*k;
mr=(700+af)/k;
if ~nargout
    plot(frq,mel,'-x');
    xlabel(['Frequency (' xticksi 'Hz)']);
    ylabel(['Frequency (' yticksi 'Mel)']);
end

function [frq,mr] = mel2frq(mel)
%MEL2FRQ  Convert Mel frequency scale to Hertz FRQ=(MEL)
%	frq = mel2frq(mel) converts a vector of Mel frequencies
%	to the corresponding real frequencies.
%   mr gives the corresponding gradients in Hz/mel.
%	The Mel scale corresponds to the perceived pitch of a tone

%	The relationship between mel and frq is given by:
%
%	m = ln(1 + f/700) * 1000 / ln(1+1000/700)
%
%  	This means that m(1000) = 1000
%
%	References:
%
%	  [1] S. S. Stevens & J. Volkman "The relation of pitch to
%		frequency", American J of Psychology, V 53, p329 1940
%	  [2] C. G. M. Fant, "Acoustic description & classification
%		of phonetic units", Ericsson Tchnics, No 1 1959
%		(reprinted in "Speech Sounds & Features", MIT Press 1973)
%	  [3] S. B. Davis & P. Mermelstein, "Comparison of parametric
%		representations for monosyllabic word recognition in
%		continuously spoken sentences", IEEE ASSP, V 28,
%		pp 357-366 Aug 1980
%	  [4] J. R. Deller Jr, J. G. Proakis, J. H. L. Hansen,
%		"Discrete-Time Processing of Speech Signals", p380,
%		Macmillan 1993
%	  [5] HTK Reference Manual p73
%	



%      Copyright (C) Mike Brookes 1998
%      Version: $Id: mel2frq.m,v 1.7 2010/08/01 08:41:57 dmb Exp $
%
%   VOICEBOX is a MATLAB toolbox for speech processing.
%   Home page: http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You can obtain a copy of the GNU General Public License from
%   http://www.gnu.org/copyleft/gpl.html or by writing to
%   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
persistent k
if isempty(k)
    k=1127.01048;
end
frq=700*sign(mel).*(exp(abs(mel)/k)-1);
mr=(700+abs(frq))/k;
if ~nargout
    plot(mel,frq,'-x');
    xlabel(['Frequency (' xticksi 'Mel)']);
    ylabel(['Frequency (' yticksi 'Hz)']);
end


function [b,c] = frq2bark(f,m)
%FRQ2BARK  Convert Hertz to BARK frequency scale BARK=(FRQ)
%       bark = frq2bark(frq) converts a vector of frequencies (in Hz)
%       to the corresponding values on the BARK scale.
% Inputs: f  matrix of frequencies in Hz
%         m  mode options
%            'h'   use high frequency correction from [1]
%            'l'   use low frequency correction from [1]
%            'H'   do not apply any high frequency correction
%            'L'   do not apply any low frequency correction
%            'z'   use the expressions from Zwicker et al. (1980) for b and c
%            's'   use the expression from Schroeder et al. (1979)
%            'u'   unipolar version: do not force b to be an odd function
%                  This has no effect on the default function which is odd anyway
%            'g'   plot a graph
%
% Outputs: b  bark values
%          c  Critical bandwidth: d(freq)/d(bark)

%   The Bark scale is named in honour of Barkhausen, the creator
%   of the unit of loudness level [2]. Criitical band k extends
%   from bark2frq(k-1) to bark2frq(k).
%
%   There are many published formulae approximating the Bark scale.
%   The default is the one from [1] but with a correction at high and
%   low frequencies to give a better fit to [2] with a continuous derivative
%   and ensure that 0 Hz = 0 Bark.
%   The h and l mode options apply the corrections from [1] which are
%   not as good and do not give a continuous derivative. The H and L
%   mode options suppress the correction entirely to give a simple formula.
%   The 's' option uses the less accurate formulae from [3] which have been
%   widely used in the lterature.
%   The 'z' option uses the formulae from [4] in which the c output
%   is not exactly the reciprocal of the derivative of the bark function.
%
%   [1] H. Traunmuller, Analytical Expressions for the
%       Tonotopic Sensory Scale�, J. Acoust. Soc. Am. 88,
%       1990, pp. 97-100.
%   [2] E. Zwicker, Subdivision of the audible frequency range into
%       critical bands, J Accoust Soc Am 33, 1961, p248.
%   [3] M. R. Schroeder, B. S. Atal, and J. L. Hall. Optimizing digital
%       speech coders by exploiting masking properties of the human ear.
%       J. Acoust Soc Amer, 66 (6): 1647�1652, 1979. doi: 10.1121/1.383662.
%   [4] E. Zwicker and E. Terhardt.  Analytical expressions for
%       critical-band rate and critical bandwidth as a function of frequency.
%       J. Acoust Soc Amer, 68 (5): 1523�1525, Nov. 1980.

%   The following code reproduces the graphs 3(c) and 3(d) from [1].
%       b0=(0:0.5:24)';
%       f0=[[2 5 10 15 20 25 30 35 40 45 51 57 63 70 77 ...
%           84 92 100 108 117 127 137 148 160 172 185 200 ...
%           215 232 250 270 290 315]*10 [34 37 40 44 48 53 ...
%           58 64 70 77 85 95 105 120 135 155]*100]';
%       b1=frq2bark(f0);      b2=frq2bark(f0,'lh');
%       b3=frq2bark(f0,'LH'); b4=frq2bark(f0,'z');
%       plot(b0,[b0 b1 b2 b3 b4]-repmat(b0,1,5));
%       xlabel('Frequency (Bark)'); ylabel('Error (Bark)');
%       legend('Exact','voicebox','Traunmuller1990', ...
%              'Traunmuller1983','Zwicker1980','Location','South');

%      Copyright (C) Mike Brookes 2006-2010
%      Version: $Id: frq2bark.m,v 1.8 2010/01/06 10:12:04 dmb Exp $
%
%   VOICEBOX is a MATLAB toolbox for speech processing.
%   Home page: http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You can obtain a copy of the GNU General Public License from
%   http://www.gnu.org/copyleft/gpl.html or by writing to
%   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
persistent A B C D P Q R S T U
if isempty(P)
    A=26.81;
    B=1960;
    C=-0.53;
    D=A*B;
    P=(0.53/(3.53)^2);
    Q=0.25;
    R=20.4;
    xy=2;
    S=0.5*Q/xy;
    T=R+0.5*xy;
    U=T-xy;
end
if nargin<2
    m=' ';
end
if any(m=='u')
    g=f;
else
    g=abs(f);
end
if any(m=='z')
    b=13*atan(0.00076*g)+3.5*atan((f/7500).^2);
    c=25+75*(1+1.4e-6*f.^2).^0.69;
elseif any(m=='s')
    b=7*log(g/650+sqrt(1+(g/650).^2));
    c=cosh(b/7)*650/7;
else
    b=A*g./(B+g)+C;
    d=D*(B+g).^(-2);
    if any(m=='l')
        m1=(b<2);
        d(m1)=d(m1)*0.85;
        b(m1)=0.3+0.85*b(m1);
    elseif ~any(m=='L')
        m1=(b<3);
        b(m1)=b(m1)+P*(3-b(m1)).^2;
        d(m1)=d(m1).*(1-2*P*(3-b(m1)));
    end
    if any(m=='h')
        m1=(b>20.1);
        d(m1)=d(m1)*1.22;
        b(m1)=1.22*b(m1)-4.422;
    elseif ~any(m=='H')
        m2=(b>T);
        m1=(b>U) & ~m2;
        b(m1)=b(m1)+S*(b(m1)-U).^2;
        b(m2)=(1+Q)*b(m2)-Q*R;
        d(m2)=d(m2).*(1+Q);
        d(m1)=d(m1).*(1+2*S*(b(m1)-U));
    end
    c=d.^(-1);
end
if ~any(m=='u')
    b=b.*sign(f);          % force to be odd
end

if ~nargout || any(m=='g')
    subplot(212)
    semilogy(f,c,'-r');
    ha=gca;
    ylabel(['Critical BW (' yticksi 'Hz)']);
    xlabel(['Frequency (' xticksi 'Hz)']);
    subplot(211)
    plot(f,b,'x-b');
    hb=gca;
    ylabel('Bark');
    xlabel(['Frequency (' xticksi 'Hz)']);
    linkaxes([ha hb],'x');
end

function [erb,bnd] = frq2erb(frq)
%FRQ2ERB  Convert Hertz to ERB frequency scale ERB=(FRQ)
%	erb = frq2erb(frq) converts a vector of frequencies (in Hz)
%	to the corresponding values on the ERB-rate scale on which
%  	the human ear has roughly constant resolution as judged by
%  	psychophysical measurements of the cochlear filters.

%	We have df/de = 6.23*f^2 + 93.39*f + 28.52
%	where the above expression gives the Equivalent Rectangular
%	Bandwidth (ERB)in Hz  of a human auditory filter with a centre
%	frequency of f kHz.
%
%	By integrating the reciprocal of the above expression, we
%	get:
%		e = a ln((f/p-1)/(f/q-1))
%
%	where p and q are the roots of the equation: -0.312 and -14.7
%  	and a = 1000/(6.23*(p-q)) = 11.17268
%
%	We actually implement e as
%
%		e = a ln (1 + b*f/(f+c))
%
%	where b = q/p - 1 = 46.06538
%	      c = -1000q = 14678.49
%	and f is in Hz
%
%	References:
%
%	  [1] B.C.J.Moore & B.R.Glasberg "Suggested formula for
%		calculating auditory-filter bandwidth and excitation
%		patterns", J Acoust Soc America V74, pp 750-753, 1983
%
%	  [2] O. Ghitza, "Auditory Models & Human Performance in Tasks
%		related to Speech Coding & Speech Recognition",
%		IEEE Trans on Speech & Audio Processing, Vol 2,
%		pp 115-132, Jan 1994
%	



%      Copyright (C) Mike Brookes 1998
%      Version: $Id: frq2erb.m,v 1.6 2010/07/18 18:39:04 dmb Exp $
%
%   VOICEBOX is a MATLAB toolbox for speech processing.
%   Home page: http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You can obtain a copy of the GNU General Public License from
%   http://www.gnu.org/copyleft/gpl.html or by writing to
%   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g=abs(frq);
erb=11.17268*sign(frq).*log(1+46.06538*g./(g+14678.49));
bnd=6.23e-6*g.^2 + 93.39e-3*g + 28.52;
if ~nargout
    plot(frq,erb,'-x');
    xlabel(['Frequency (' xticksi 'Hz)']);
    ylabel(['Frequency (' yticksi 'Erb-rate)']);
end


function [f,t]=enframe(x,win,inc)
%ENFRAME split signal up into (overlapping) frames: one per row. [F,T]=(X,WIN,INC)
%
%	F = ENFRAME(X,LEN) splits the vector X(:) up into
%	frames. Each frame is of length LEN and occupies
%	one row of the output matrix. The last few frames of X
%	will be ignored if its length is not divisible by LEN.
%	It is an error if X is shorter than LEN.
%
%	F = ENFRAME(X,LEN,INC) has frames beginning at increments of INC
%	The centre of frame I is X((I-1)*INC+(LEN+1)/2) for I=1,2,...
%	The number of frames is fix((length(X)-LEN+INC)/INC)
%
%	F = ENFRAME(X,WINDOW) or ENFRAME(X,WINDOW,INC) multiplies
%	each frame by WINDOW(:)
%
%   The second output argument, T, gives the time in samples at the centre
%   of each frame. T=i corresponds to the time of sample X(i). 
%
% Example of frame-based processing:
%          INC=20       													% set frame increment
%          NW=INC*2     													% oversample by a factor of 2 (4 is also often used)
%          S=cos((0:NW*7)*6*pi/NW);								% example input signal
%          W=sqrt(hamming(NW+1)); W(end)=[];      % sqrt hamming window of period NW
%          F=enframe(S,W,INC);               			% split into frames
%          ... process frames ...
%          X=overlapadd(F,W,INC);           			% reconstitute the time waveform (omit "X=" to plot waveform)

%	   Copyright (C) Mike Brookes 1997
%      Version: $Id: enframe.m,v 1.7 2009/11/01 21:08:21 dmb Exp $
%
%   VOICEBOX is a MATLAB toolbox for speech processing.
%   Home page: http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You can obtain a copy of the GNU General Public License from
%   http://www.gnu.org/copyleft/gpl.html or by writing to
%   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nx=length(x(:));
nwin=length(win);
if (nwin == 1)
   len = win;
else
   len = nwin;
end
if (nargin < 3)
   inc = len;
end
nf = fix((nx-len+inc)/inc);
f=zeros(nf,len);
indf= inc*(0:(nf-1)).';
inds = (1:len);
f(:) = x(indf(:,ones(1,len))+inds(ones(nf,1),:));
if (nwin > 1)
    w = win(:)';
    f = f .* w(ones(nf,1),:);
end
if nargout>1
    t=(1+len)/2+indf;
end

function y=rfft(x,n,d)
%RFFT     Calculate the DFT of real data Y=(X,N,D)
% Data is truncated/padded to length N if specified.
%   N even:	(N+2)/2 points are returned with
% 			the first and last being real
%   N odd:	(N+1)/2 points are returned with the
% 			first being real
% In all cases fix(1+N/2) points are returned
% D is the dimension along which to do the DFT



%      Copyright (C) Mike Brookes 1998
%      Version: $Id: rfft.m,v 1.8 2010/12/10 14:25:00 dmb Exp $
%
%   VOICEBOX is a MATLAB toolbox for speech processing.
%   Home page: http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You can obtain a copy of the GNU General Public License from
%   http://www.gnu.org/copyleft/gpl.html or by writing to
%   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

s=size(x);
if prod(s)==1
    y=x
else
    if nargin <3 || isempty(d)
        d=find(s>1,1);
        if nargin<2
            n=s(d);
        end
    end
    if isempty(n) 
        n=s(d);
    end
    y=fft(x,n,d);
    y=reshape(y,prod(s(1:d-1)),n,prod(s(d+1:end))); 
    s(d)=1+fix(n/2);
    y(:,s(d)+1:end,:)=[];
    y=reshape(y,s);
end

function [x,cf,il,ih]=filtbankm(p,n,fs,fl,fh,w)
%FILTBANKM determine matrix for a linear/mel/erb/bark-spaced filterbank [X,MN,MX]=(P,N,FS,FL,FH,W)
%
% Inputs:
%       p   number of filters in filterbank or the filter spacing in k-mel/bark/erb [ceil(4.6*log10(fs))]
%		n   length of fft
%           or [nfrq dfrq frq1] nfrq=number of input frequency bins, frequency increment (Hz), first bin freq (Hz)
%		fs  sample rate in Hz
%		fl  low end of the lowest filter in Hz (see 'h' option) [default = 0 or 30Hz for 'l' option]
%		fh  high end of highest filter in Hz [default = fs/2]
%		w   any sensible combination of the following:
%
%             'b' = bark scale instead of mel
%             'e' = erb-rate scale
%             'l' = log10 Hz frequency scale
%             'f' = linear frequency scale [default]
%             'm' = mel frequency scale
%
%             'c' = fl & fh specify centre of low and high filters instead of edges
%             'h' = fl & fh are in mel/erb/bark/log10 instead of Hz
%             'H' = cf outputs are in mel/erb/bark/log10 instead of Hz
%
%		      'y' = lowest filter remains at 1 down to 0 frequency and
%			        highest filter remains at 1 up to nyquist freqency
%		            The total power in the fft is preserved (unless 'u' is specified).
%             'Y' = extend only at low frequency end (or high end if 'y' also specified)
%
%             'p' = input P specifies the number of filters [default if P>=1]
%             'P' = input P specifies the filter spacing [default if P<1]
%
%             'u' = input and output are power per Hz instead of power.
%             'U' = input is power but output is power per Hz.
%
%             's' = single-sided input: do not include symmetric negative frequencies (i.e. non-DC inputs have been doubled)
%             'S' = single-sided output: do not mirror the non-DC filter characteristics (i.e. double non-DC outputs)
%
%             'g' = plot filter coefficients as graph
%             'G' = plot filter coefficients as image [default if no output arguments present]
%
%
% Outputs:	x     a sparse matrix containing the filterbank amplitudes
%		          If the il and ih outputs are given then size(x)=[p,mx-mn+1]
%                 otherwise size(x)=[p,1+floor(n/2)]
%                 Note that the peak filter values equal 2 to account for the power
%                 in the negative FFT frequencies.
%           cf    the filterbank centre frequencies in Hz (see 'H' option)
%		    il    the lowest fft bin with a non-zero coefficient
%		    ih    the highest fft bin with a non-zero coefficient
%
% The routine performs interpolation of the input spectrum by convolving the power spectrum
% with a triangular filter and then simulates a filterbank with asymetric triangular filters.
%
% Examples of use:
%
% (a) Calcuate the Mel-frequency Cepstral Coefficients
%
%       f=rfft(s);			        % rfft() returns only 1+floor(n/2) coefficients
%		x=filtbankm(p,n,fs,0,fs/2,'m');	        % n is the fft length, p is the number of filters wanted
%		z=log(x*abs(f).^2);         % multiply x by the power spectrum
%		c=dct(z);                   % take the DCT
%
% (b) Calcuate the Mel-frequency Cepstral Coefficients efficiently
%
%       f=fft(s);                        % n is the fft length, p is the number of filters wanted
%       [x,cf,na,nb]=filtbankm(p,n,fs,0,fs/2,'m');   % na:nb gives the fft bins that are needed
%       z=log(x*(f(na:nb)).*conj(f(na:nb)));
%
% (c) Plot the calculated filterbanks
%
%      plot((0:floor(n/2))*fs/n,filtbankm(p,n,fs,0,fs/2,'m')')   % fs=sample frequency
%
% (d) Plot the filterbanks
%
%      filtbankm(p,n,fs,0,fs/2,'m');
%
% References:
%
% [1] S. S. Stevens, J. Volkman, and E. B. Newman. A scale for the measurement
%     of the psychological magnitude of pitch. J. Acoust Soc Amer, 8: 185�19, 1937.
% [2] S. Davis and P. Mermelstein. Comparison of parametric representations for
%     monosyllabic word recognition in continuously spoken sentences.
%     IEEE Trans Acoustics Speech and Signal Processing, 28 (4): 357�366, Aug. 1980.

% Bugs/Suggestions
% (1) default frequencies won't work if the h option is specified
% (2) low default frequency is invalid if the 'l' option is specified

%      Copyright (C) Mike Brookes 1997-2009
%      Version: $Id: filtbankm.m,v 1.5 2011/04/20 16:08:35 dmb Exp $
%
%   VOICEBOX is a MATLAB toolbox for speech processing.
%   Home page: http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You can obtain a copy of the GNU General Public License from
%   http://www.gnu.org/copyleft/gpl.html or by writing to
%   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Note "FFT bin_0" assumes DC = bin 0 whereas "FFT bin_1" means DC = bin 1

if nargin < 6
    w='f'; % default option: linear frequency scale
end
wr=' ';   % default warping is linear frequency
for i=1:length(w)
    if any(w(i)=='lebm');
        wr=w(i);
    end
end
if nargin < 5 || ~numel(fh)
    fh=0.5*fs; % max freq is the nyquist
end
if nargin < 4 || ~numel(fl)
    if wr=='l'
        fl=30;  % min freq is 30 Hz for log scale
    else
    fl=0; % min freq is DC
    end
end

f1=0;
if numel(n)>1
    nf=n(1);  % number of input frequency bins
    df=n(2);  % input frequency bin spacing
    if numel(n)>2
        f1=n(3); % frequency of first bin
    end
else
    nf=1+floor(n/2); % number of input frequency bins
    df=fs/n;  % input frequency bin spacing
end
fin0=f1+(0:nf-1)*df;  % input frequency bins

mflh=[fl fh];
if ~any(w=='h')             % convert Hz to mel/erb/...
    switch wr
        case 'm'
            mflh=frq2mel(mflh);       % convert frequency limits into mel
        case 'l'
            if fl<=0
                error('Low frequency limit must be >0 for l option');
            end
            mflh=log10(mflh);       % convert frequency limits into log10 Hz
        case 'e'
            mflh=frq2erb(mflh);       % convert frequency limits into erb-rate
        case 'b'
            mflh=frq2bark(mflh);       % convert frequency limits into bark
    end
end
melrng=mflh*(-1:2:1)';          % mel/erb/... range
% fn2=floor(n/2);     % bin index of highest positive frequency (Nyquist if n is even)
if isempty(p)
    p=ceil(4.6*log10(2*(f1+(nf-1)*df)));         % default number of output filters
end
puc=any(w=='P') || (p<1) && ~any(w=='p');
if any(w=='c')              % c option: specify fiter centres not edges
    if puc
        p=round(melrng/(p*1000))+1;
    end
    melinc=melrng/(p-1);
    mflh=mflh+(-1:2:1)*melinc;
else
    if puc
        p=round(melrng/(p*1000))-1;
    end
    melinc=melrng/(p+1);
end
%
% Calculate the FFT bins0 corresponding to the filters
%
cf=mflh(1)+(0:p+1)*melinc; % centre frequencies in mel/erb/... including dummy ends
cf(2:end)=max(cf(2:end),0); % only the first point can be negative
switch wr    % convert centre frequencies from mel/erb/... to Hz
    case 'l'
        mb=10.^(cf);
    case 'e'
        mb=erb2frq(cf);
    case 'b'
        mb=bark2frq(cf);
    case 'm'
        mb=mel2frq(cf);
    otherwise
        mb=cf;
end

% first sort out 2-sided input frequencies

fin=fin0;
fin(nf+1)=fin(nf)+df; % add on a dummy point at the high end
if fin(1)==0
    fin=[-fin(nf+1:-1:2) fin];
elseif fin(1)<=df/2
    fin=[-fin(nf+1:-1:1) fin];
elseif fin(1)<df
    fin=[-fin(nf+1:-1:1) fin(1)-df df-fin(1) fin];
elseif fin(1)==df
    fin=[-fin(nf+1:-1:1) 0 fin];
else
    fin=[-fin(nf+1:-1:1) df-fin(1) fin(1)-df fin];
end
nfin=length(fin);  % length of extended input frequency list

% now sort out the interleaving

fout=mb;  % output frequencies in Hz
lowex=any(w=='y')~=any(w=='Y');   % extend to 0 Hz
highex=any(w=='y') && (fout(end-1)<fin(end));  % extend at high end
if lowex
    fout=[0 0 fout(2:end)];
end
if highex
    fout=[fout(1:end-1) fin(end) fin(end)];
end
mfout=length(fout);
if any(w=='u') || any(w=='U')
    gout=fout(3:mfout)-fout(1:mfout-2);
    gout=2*(gout+(gout==0)).^(-1); % Gain of output triangles
else
    gout=ones(1,mfout-2);
end
if any(w=='S')
    msk=fout(2:mfout-1)~=0;
    gout(msk)=2*gout(msk); % double non-DC outputs for a 1-sided ouptu spectrum
end
if any(w=='u')
    gin=ones(1,nfin-2);
else
    gin=fin(3:nfin)-fin(1:nfin-2);
    gin=2*(gin+(gin==0)).^(-1); % Gain of input triangles
end
msk=fin(2:end-1)==0;
if any(w=='s')
    gin(~msk)=0.5*gin(~msk); % halve non-DC inputs to change back to a 2-sided spectrum
end
if lowex
    gin(msk)=2*gin(msk);  % double DC input to preserve its power
end
foutin=[fout fin];
nfall=length(foutin);
wleft=[0 fout(2:mfout)-fout(1:mfout-1) 0 fin(2:nfin)-fin(1:nfin-1)]; % left width
wright=[wleft(2:end) 0]; % right width
ffact=[0 gout 0 0 gin(1:min(nf,nfin-nf-2)) zeros(1,max(nfin-2*nf-2,0)) gin(nfin-nf-1:nfin-2) 0]; % gain of triangle posts
% ffact(wleft+wright==0)=0; % disable null width triangles shouldn't need this if all frequencies are distinct
[fall,ifall]=sort(foutin);
jfall=zeros(1,nfall);
infall=1:nfall;
jfall(ifall)=infall; % unsort->sort index
ffact(ifall([1:max(jfall(1),jfall(mfout+1))-2 min(jfall(mfout),jfall(nfall))+2:nfall]))=0;  % zap nodes that are much too small/big

nxto=cumsum(ifall<=mfout);
nxti=cumsum(ifall>mfout);
nxtr=min(nxti+1+mfout,nfall);  % next input node to the right of each value (or nfall if none)
nxtr(ifall>mfout)=1+nxto(ifall>mfout); % next post to the right of opposite type (unsorted indexes)
nxtr=nxtr(jfall);  % next post to the right of opposite type (converted to unsorted indices) or if none: nfall/(mfout+1)

% each triangle is "attached" to the node at its extreme right end
% the general result for integrating the product of two trapesiums with
% heights (a,b) and (c,d) over a width x is (ad+bc+2bd+2ac)*w/6
%
% integrate product of lower triangles

msk0=(ffact>0);
msk=msk0 & (ffact(nxtr)>0); % select appropriate triangle pairs (unsorted indices)
ix1=infall(msk); % unsorted indices of leftmost post of pair
jx1=nxtr(msk);  % unsorted indices of rightmost post of pair
vfgx=foutin(ix1)-foutin(jx1-1); % length of right triangle to the left of the left post
yx=min(wleft(ix1),vfgx); % integration length
wx1=ffact(ix1).*ffact(jx1).*yx.*(wleft(ix1).*vfgx-yx.*(0.5*(wleft(ix1)+vfgx)-yx/3))./(wleft(ix1).*wleft(jx1)+(yx==0));

% integrate product of upper triangles

nxtu=max([nxtr(2:end)-1 0],1);
msk=msk0 & (ffact(nxtu)>0);
ix2=infall(msk); % unsorted indices of leftmost post of pair
jx2=nxtu(msk);  % unsorted indices of rightmost post of pair
vfgx=foutin(ix2+1)-foutin(jx2); % length of left triangle to the right of the right post
yx=min(wright(ix2),vfgx); % integration length
yx(foutin(jx2+1)<foutin(ix2+1))=0; % zap invalid triangles
wx2=ffact(ix2).*ffact(jx2).*yx.^2.*((0.5*(wright(jx2)-vfgx)+yx/3))./(wright(ix2).*wright(jx2)+(yx==0));

% integrate lower triangle and upper triangle that ends to its right

nxtu=max(nxtr-1,1);
msk=msk0 & (ffact(nxtu)>0);
ix3=infall(msk); % unsorted indices of leftmost post of pair
jx3=nxtu(msk);  % unsorted indices of rightmost post of pair
vfgx=foutin(ix3)-foutin(jx3); % length of upper triangle to the left of the lower post
yx=min(wleft(ix3),vfgx); % integration length
yx(foutin(jx3+1)<foutin(ix3))=0; % zap invalid triangles
wx3=ffact(ix3).*ffact(jx3).*yx.*(wleft(ix3).*(wright(jx3)-vfgx)+yx.*(0.5*(wleft(ix3)-wright(jx3)+vfgx)-yx/3))./(wleft(ix3).*wright(jx3)+(yx==0));

% integrate upper triangle and lower triangle that starts to its right

nxtu=[nxtr(2:end) 1];
msk=msk0 & (ffact(nxtu)>0);
ix4=infall(msk); % unsorted indices of leftmost post of pair
jx4=nxtu(msk);  % unsorted indices of rightmost post of pair
vfgx=foutin(ix4+1)-foutin(jx4-1); % length of upper triangle to the left of the lower post
yx=min(wright(ix4),vfgx); % integration length
wx4=ffact(ix4).*ffact(jx4).*yx.^2.*(0.5*vfgx-yx/3)./(wright(ix4).*wleft(jx4)+(yx==0));

% now create the matrix

iox=sort([ix1 ix2 ix3 ix4;jx1 jx2 jx3 jx4]);
% [iox;[wx1 wx2 wx3 wx4]>0 ]
msk=iox(2,:)<=(nfall+mfout)/2;
iox(2,msk)=(nfall+mfout+1)-iox(2,msk);  % convert negative frequencies to positive
if highex
    iox(1,iox(1,:)==mfout-1)=mfout-2; % merge highest two output nodes
end
if lowex
    iox(1,iox(1,:)==2)=3; % merge lowest two output nodes
end

x=sparse(iox(1,:)-1-lowex,max(iox(2,:)-nfall+nf+1,1),[wx1 wx2 wx3 wx4],p,nf);
%
% sort out the output argument options
%
if ~any(w=='H')
    cf=mb;         % output Hz instead of mel/erb/...
end
cf=cf(2:p+1);  % remove dummy end frequencies
il=1;
ih=nf;
if nargout > 2
    msk=full(any(x>0,1));
    il=find(msk,1);
    if ~numel(il)
        ih=1;
    elseif nargout >3
        ih=find(msk,1,'last');
    end
    x=x(:,il:ih);
end
if any(w=='u')
    sx=sum(x,2);
    x=x./repmat(sx+(sx==0),1,size(x,2));
end
%
% plot results if no output arguments or g option given
%
if ~nargout || any(w=='g') || any(w=='G') % plot idealized filters
    if ~any(w=='g') && ~any(w=='G')
        w=[w 'G'];
    end
    newfig=0;
    if  any(w=='g')
        plot(f1-df+(il:ih)*df,x');
        title(['filtbankm: mode = ' w]);
        xlabel(['Frequency (' xticksi 'Hz)']);
        ylabel('Weight');
        newfig=1;
    end

    if  any(w=='G')
        if newfig
            figure;
        end
        imagesc(f1-df+(il:ih)*df,1:p,x);
        axis 'xy'
        colorbar;
        cblabel('Weight');
        switch wr
            case 'l'
                type='Log-spaced';
            case 'e'
                type='Erb-spaced';
            case 'b'
                type='Bark-spaced';
            case 'm'
                type='Mel-spaced';
            otherwise
                type='Linear-spaced';
        end
        ylabel([type ' Filter']);
        xlabel(['Frequency (' xticksi 'Hz)']);
        title(['filtbankm: mode = ' w]);
    end

end

function [frq,bnd] = erb2frq(erb)
%ERB2FRQ  Convert ERB frequency scale to Hertz FRQ=(ERB)
%	frq = erb2frq(erb) converts a vector of ERB-rate values
%	to the corresponding frequencies in Hz.
%   [frq,bnd] =  erb2frq(erb) also calculates the ERB bandwidths
%
%       Note that erb values must not exceed 42.79
%
% See also: frq2erb

%	We have df/de = 6.23*f^2 + 93.39*f + 28.52
%	where the above expression gives the Equivalent Rectangular
%	Bandwidth (ERB)in Hz  of a human auditory filter with a centre
%	frequency of f kHz.
%
%	By integrating the reciprocal of the above expression, we
%	get:
%		e = a ln((f/p-1)/(f/q-1))/c
%
%	where p and q are the roots of the equation: -0.312 and -14.7
%  	and c = (6.23*(p-q))/1000 = 0.08950404
%
%	from this we can derive:
%
%	f = a/(b-exp(c*e)) - d
%
%	where a = 1000 q (1 - q/p) = 676170.4
%	      b = q/p = 47.06538
%	      d = -1000q = 14678.49
%	and f is in Hz
%
%	References:
%
%	  [1] B.C.J.Moore & B.R.Glasberg "Suggested formula for
%		calculating auditory-filter bandwidth and excitation
%		patterns", J Acoust Soc America V74, pp 750-753, 1983
%
%	  [2] O. Ghitza, "Auditory Models & Human Performance in Tasks
%		related to Speech Coding & Speech Recognition",
%		IEEE Trans on Speech & Audio Processing, Vol 2,
%		pp 115-132, Jan 1994
%	

%      Copyright (C) Mike Brookes 1998
%      Version: $Id: erb2frq.m 713 2011-10-16 14:45:43Z dmb $
%
%   VOICEBOX is a MATLAB toolbox for speech processing.
%   Home page: http://www.ee.ic.ac.uk/hp/staff/dmb/voicebox/voicebox.html
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   This program is free software; you can redistribute it and/or modify
%   it under the terms of the GNU General Public License as published by
%   the Free Software Foundation; either version 2 of the License, or
%   (at your option) any later version.
%
%   This program is distributed in the hope that it will be useful,
%   but WITHOUT ANY WARRANTY; without even the implied warranty of
%   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%   GNU General Public License for more details.
%
%   You can obtain a copy of the GNU General Public License from
%   http://www.gnu.org/copyleft/gpl.html or by writing to
%   Free Software Foundation, Inc.,675 Mass Ave, Cambridge, MA 02139, USA.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

frq = sign(erb).*(676170.4*(47.06538-exp(0.08950404*abs(erb))).^(-1) - 14678.49);
bnd=6.23e-6*frq.^2 + 93.39e-3*abs(frq) + 28.52;
if ~nargout
    plot(erb,frq,'-x');
    xlabel(['Frequency (' xticksi 'Erb-rate)']);
    ylabel(['Frequency (' yticksi 'Hz)']);
end
