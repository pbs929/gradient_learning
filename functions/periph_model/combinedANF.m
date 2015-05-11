function [gamma, f, ANF] =combinedANF(sw, fs, numChannels, highFreq, lowFreq)
% This function combines the various elements of the Slaney package with
% parameters optimized for processing birdsongs

% --load sounds--
if size(sw,2)>size(sw,1), sw=sw'; end
normfactor=0.5*10^4;
sw = sw*normfactor;

% --Gammatone filter--
[~,~,f]=MakeERBFilters_printB(2*highFreq,numChannels,lowFreq);
[fcoefs,b,f]=MakeERBFilters_printB(fs,f,lowFreq);
b=flipud(b); f=flipud(f); fcoefs=flipud(fcoefs);

gamma = ERBFilterBank(sw, fcoefs);

if nargout < 3, return; end

% --Meddis Hair Cell-- - parameters from a later paper
M = 1;
A = 5;
B = 800;
g = 1000;
y = 5.05;
l = 1250;
r = 6580;
x = 66.31;
h = 50000;

ANF=MeddisHairCell_params(gamma,fs,M,A,B,g,y,l,r,x,h);

end

function [fcoefs,b,cf]=MakeERBFilters_printB(fs,numChannels,lowFreq)
% function [fcoefs]=MakeERBFilters(fs,numChannels,lowFreq)
% This function computes the filter coefficients for a bank of 
% Gammatone filters.  These filters were defined by Patterson and 
% Holdworth for simulating the cochlea.  
% 
% The result is returned as an array of filter coefficients.  Each row 
% of the filter arrays contains the coefficients for four second order 
% filters.  The transfer function for these four filters share the same
% denominator (poles) but have different numerators (zeros).  All of these
% coefficients are assembled into one vector that the ERBFilterBank 
% can take apart to implement the filter.
%
% The filter bank contains "numChannels" channels that extend from
% half the sampling rate (fs) to "lowFreq".  Alternatively, if the numChannels
% input argument is a vector, then the values of this vector are taken to
% be the center frequency of each desired filter.  (The lowFreq argument is
% ignored in this case.)

% Note this implementation fixes a problem in the original code by
% computing four separate second order filters.  This avoids a big
% problem with round off errors in cases of very small cfs (100Hz) and
% large sample rates (44kHz).  The problem is caused by roundoff error
% when a number of poles are combined, all very close to the unit
% circle.  Small errors in the eigth order coefficient, are multiplied
% when the eigth root is taken to give the pole location.  These small
% errors lead to poles outside the unit circle and instability.  Thanks
% to Julius Smith for leading me to the proper explanation.

% Execute the following code to evaluate the frequency
% response of a 10 channel filterbank.
%	fcoefs = MakeERBFilters(16000,10,100);
%	y = ERBFilterBank([1 zeros(1,511)], fcoefs);
%	resp = 20*log10(abs(fft(y')));
%	freqScale = (0:511)/512*16000;
%	semilogx(freqScale(1:255),resp(1:255,:));
%	axis([100 16000 -60 0])
%	xlabel('Frequency (Hz)'); ylabel('Filter Response (dB)');

% Rewritten by Malcolm Slaney@Interval.  June 11, 1998.
% (c) 1998 Interval Research Corporation  

T = 1/fs;
if length(numChannels) == 1
	cf = ERBSpace(lowFreq, fs/2, numChannels);
else
	cf = numChannels(1:end);
	if size(cf,2) > size(cf,1)
		cf = cf';
	end
end

% Change the following three parameters if you wish to use a different
% ERB scale.  Must change in ERBSpace too.
EarQ = 9.26449;				%  Glasberg and Moore Parameters
minBW = 24.7;
order = 1;

ERB = ((cf/EarQ).^order + minBW^order).^(1/order);  %64vector (row vector) - ERB for each band
B=1.019*2*pi*ERB;  %64vector - b values for each band (assumes 4th order) B=2*pi*b (b=1.019ERB)
b=1.019*ERB;   % MY ADDITION - display the value of b for the gammatone.  

A0 = T;  %scalar - sampling period
A2 = 0;  %scalar - zero
B0 = 1;  %scalar - one
B1 = -2*cos(2*cf*pi*T)./exp(B*T);  %64vector - resembles part of gammatone?
B2 = exp(-2*B*T);  %64vector

A11 = -(2*T*cos(2*cf*pi*T)./exp(B*T) + 2*sqrt(3+2^1.5)*T*sin(2*cf*pi*T)./ ...
		exp(B*T))/2;  %64vector
A12 = -(2*T*cos(2*cf*pi*T)./exp(B*T) - 2*sqrt(3+2^1.5)*T*sin(2*cf*pi*T)./ ...
		exp(B*T))/2;  %64vector
A13 = -(2*T*cos(2*cf*pi*T)./exp(B*T) + 2*sqrt(3-2^1.5)*T*sin(2*cf*pi*T)./ ...
		exp(B*T))/2;  %64vector
A14 = -(2*T*cos(2*cf*pi*T)./exp(B*T) - 2*sqrt(3-2^1.5)*T*sin(2*cf*pi*T)./ ...
		exp(B*T))/2;  %64vector

gain = abs((-2*exp(4*i*cf*pi*T)*T + ...
                 2*exp(-(B*T) + 2*i*cf*pi*T).*T.* ...
                         (cos(2*cf*pi*T) - sqrt(3 - 2^(3/2))* ...
                          sin(2*cf*pi*T))) .* ...
           (-2*exp(4*i*cf*pi*T)*T + ...
             2*exp(-(B*T) + 2*i*cf*pi*T).*T.* ...
              (cos(2*cf*pi*T) + sqrt(3 - 2^(3/2)) * ...
               sin(2*cf*pi*T))).* ...
           (-2*exp(4*i*cf*pi*T)*T + ...
             2*exp(-(B*T) + 2*i*cf*pi*T).*T.* ...
              (cos(2*cf*pi*T) - ...
               sqrt(3 + 2^(3/2))*sin(2*cf*pi*T))) .* ...
           (-2*exp(4*i*cf*pi*T)*T + 2*exp(-(B*T) + 2*i*cf*pi*T).*T.* ...
           (cos(2*cf*pi*T) + sqrt(3 + 2^(3/2))*sin(2*cf*pi*T))) ./ ...
          (-2 ./ exp(2*B*T) - 2*exp(4*i*cf*pi*T) +  ...
           2*(1 + exp(4*i*cf*pi*T))./exp(B*T)).^4);  %64vector
	
allfilts = ones(length(cf),1);  %64vector of ones
fcoefs = [A0*allfilts A11 A12 A13 A14 A2*allfilts B0*allfilts B1 B2 gain];  %array of 64vectors

if (0)						% Test Code
	A0  = fcoefs(:,1);
	A11 = fcoefs(:,2);
	A12 = fcoefs(:,3);
	A13 = fcoefs(:,4);
	A14 = fcoefs(:,5);
	A2  = fcoefs(:,6);
	B0  = fcoefs(:,7);
	B1  = fcoefs(:,8);
	B2  = fcoefs(:,9);
	gain= fcoefs(:,10);	
	chan=1;
	x = [1 zeros(1, 511)];
	y1=filter([A0(chan)/gain(chan) A11(chan)/gain(chan) ...
		A2(chan)/gain(chan)],[B0(chan) B1(chan) B2(chan)], x);
	y2=filter([A0(chan) A12(chan) A2(chan)], ...
			[B0(chan) B1(chan) B2(chan)], y1);
	y3=filter([A0(chan) A13(chan) A2(chan)], ...
			[B0(chan) B1(chan) B2(chan)], y2);
	y4=filter([A0(chan) A14(chan) A2(chan)], ...
			[B0(chan) B1(chan) B2(chan)], y3);
	semilogx((0:(length(x)-1))*(fs/length(x)),20*log10(abs(fft(y4))));
end
end

function cfArray = ERBSpace(lowFreq, highFreq, N)
% 
% This function computes an array of N frequencies uniformly spaced between
% highFreq and lowFreq on an ERB scale.  N is set to 100 if not specified.
%
% See also linspace, logspace, MakeERBCoeffs, MakeERBFilters.
%
% For a definition of ERB, see Moore, B. C. J., and Glasberg, B. R. (1983).
% "Suggested formulae for calculating auditory-filter bandwidths and
% excitation patterns," J. Acoust. Soc. Am. 74, 750-753.

if nargin < 1
	lowFreq = 100;
end

if nargin < 2
	highFreq = 44100/4;
end

if nargin < 3
	N = 100;
end

% Change the following three parameters if you wish to use a different
% ERB scale.  Must change in MakeERBCoeffs too.
EarQ = 9.26449;				%  Glasberg and Moore Parameters
minBW = 24.7;
order = 1;

% All of the followFreqing expressions are derived in Apple TR #35, "An
% Efficient Implementation of the Patterson-Holdsworth Cochlear
% Filter Bank."  See pages 33-34.
cfArray = -(EarQ*minBW) + exp((1:N)'*(-log(highFreq + EarQ*minBW) + ...
		log(lowFreq + EarQ*minBW))/N) * (highFreq + EarQ*minBW);

end

function output = ERBFilterBank(x, fcoefs)
% function output = ERBFilterBank(x, fcoefs)
% Process an input waveform with a gammatone filter bank. This function 
% takes a single sound vector, and returns an array of filter outputs, one 
% channel per row.
%
% The fcoefs parameter, which completely specifies the Gammatone filterbank,
% should be designed with the MakeERBFilters function.  If it is omitted,
% the filter coefficients are computed for you assuming a 22050Hz sampling
% rate and 64 filters regularly spaced on an ERB scale from fs/2 down to 100Hz.
%

% Malcolm Slaney @ Interval, June 11, 1998.
% (c) 1998 Interval Research Corporation  
% Thanks to Alain de Cheveigne' for his suggestions and improvements.

if nargin < 1
	error('Syntax: output_array = ERBFilterBank(input_vector[, fcoefs]);');
end

if nargin < 2
	fcoefs = MakeERBFilters(22050,64,100);
end

if size(fcoefs,2) ~= 10
	error('fcoefs parameter passed to ERBFilterBank is the wrong size.');
end

if size(x,2) < size(x,1)
	x = x';
end

A0  = fcoefs(:,1);
A11 = fcoefs(:,2);
A12 = fcoefs(:,3);
A13 = fcoefs(:,4);
A14 = fcoefs(:,5);
A2  = fcoefs(:,6);
B0  = fcoefs(:,7);
B1  = fcoefs(:,8);
B2  = fcoefs(:,9);
gain= fcoefs(:,10);	

output = zeros(size(gain,1), length(x));
for chan = 1: size(gain,1)
	y1=filter([A0(chan)/gain(chan) A11(chan)/gain(chan) ...
		   A2(chan)/gain(chan)], ...
				[B0(chan) B1(chan) B2(chan)], x);
	y2=filter([A0(chan) A12(chan) A2(chan)], ...
				[B0(chan) B1(chan) B2(chan)], y1);
	y3=filter([A0(chan) A13(chan) A2(chan)], ...
				[B0(chan) B1(chan) B2(chan)], y2);
	y4=filter([A0(chan) A14(chan) A2(chan)], ...
				[B0(chan) B1(chan) B2(chan)], y3);
	output(chan, :) = y4;
end

if 0
	semilogx((0:(length(x)-1))*(fs/length(x)),20*log10(abs(fft(output))));
end
end

function y = MeddisHairCell_params(data,sampleRate,M,A,B,g,y,l,r,x,h)
% y = MeddisHairCell(data,sampleRate)
% This function calculates Ray Meddis' hair cell model for a
% number of channels.  Data is arrayed as one channel per row.
% All channels are done in parallel (but each time step is
% sequential) so it will be much more efficient to process lots
% of channels at once.

% (c) 1998 Interval Research Corporation  
% Changed h and added comment at suggestion of Alain de Cheveigne'   12/11/98

subtractSpont=0; 

% Internal constants
dt = 1/sampleRate;
gdt = g*dt;
ydt = y*dt;
ldt = l*dt;
rdt = r*dt;
xdt = x*dt;
[numChannels dataLength] = size(data);

% Initial values 
kt = g*A/(A+B);
spont = M*y*kt/(l*kt+y*(l+r));
c = spont * ones(numChannels,1);
q = c*(l+r)/kt;
w = c*r/x;
zeroVector = zeros(numChannels,1);

% Now iterate through each time slice of the data.  Use the
% max function to implement the "if (0>" test.
y = zeros(numChannels, dataLength);
for i = 1:dataLength
 limitedSt = max(data(:,i)+A,0);
 kt = gdt*limitedSt./(limitedSt+B);
 replenish = max(ydt*(M-q),zeroVector);
 eject = kt.*q;
 loss = ldt.*c;
 reuptake = rdt.*c;
 reprocess = xdt.*w;

 q = q + replenish - eject + reprocess;
 c = c + eject - loss - reuptake;
 w = w + reuptake - reprocess;
 y(:,i) = c;
end

y = h .* y;

if (subtractSpont > 0)
 y=max(0,y-spont);
end
end