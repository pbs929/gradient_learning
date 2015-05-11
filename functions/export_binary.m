function nsamp = export_binary(whichFiles, param, saveName)
% This function loads peripheral model data using the function 
% 'load_data.m', concatenates in a large matrix, and exports to a binary 
% file for use with my c++ programs.
%
% Currently I have commented out the code for supplying a header to the 
% output file.  

% get list of speech files
files = rdir([param.sampleDir, '*.mat']);
if ~isempty(whichFiles), files = files(whichFiles); end

% display target file location
disp('saving to:')
disp(saveName)

% open target file
fid = fopen(saveName,'wb');

% write number of sample vectors as a placeholder
%nsamp = 40;
%fwrite(fid, nsamp, 'integer*4',0,'b');
%nsamp = 0;

% write sample size (dimension of vedctor)
%ndim = 256;
%fwrite(fid,ndim,'integer*4',0,'b');

% loop through files
for fle = whichFiles % male files start at 56
    % load a file
    [~,~,~,data]=load_data(files(fle).name,param.sampleDir,param.periph,...
        param.normType,param.int,param.delay);
    pause
    data = data(:,param.delay:end);
    disp(files(fle).name)
    
    %[dimension,nsamples] = size(data);
    
    %nsamp = nsamp+nsamples;
    fwrite(fid,data,'real*8',0,'n');
end

%return to beginning of file and write true nsamp
%fseek(fid,0,'bof');
%fwrite(fid,nsamp,'integer*4',0,'b');

fclose(fid);
