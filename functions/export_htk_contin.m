function export_htk_contin(strfDir,param,dirname)
% This version for exporting the original continuous data samples

% make storage directory if it doesn't exist
if ~exist(dirname,'dir')
    mkdir(dirname);
end

% load the strfs
strfs = load(strfDir);

% training data
if param.runClean
    sampledir = [param.trainDir,'train_clean']; % where the original samples are
    subdir = 'train_clean/'; % subdirectory to put the outputs
    writeseqs(strfs,param,dirname,sampledir,subdir);
end

% testing data
if param.runClean
    sampledir = [param.testDir,'testa_clean']; % where the original samples are
    subdir = 'testa_clean/'; % subdirectory to put the outputs
    writeseqs(strfs,param,dirname,sampledir,subdir);
end

% noisy testing
for iLev = 1:length(param.noiseLevs)
    lev = param.noiseLevs(iLev);
    for nse = param.noiseType
        % load test file
        sampledir = [param.testDir,'n',num2str(nse),'/testa_n',num2str(nse),'_snr',num2str(lev),]; % where the original samples are
        subdir = ['testa/n',num2str(nse),'_snr',num2str(lev),'/']; % subdirectory to put the outputs
        writeseqs(strfs,param,dirname,sampledir,subdir);
    end
end

% noisy training data
sampledir = [param.trainDir,'train_multi/train']; % where the original samples are
subdir = 'train_multi/'; % subdirectory to put the outputs
writeseqs(strfs,param,dirname,sampledir,subdir);

end


function writeseqs(strf, param, dirname, sampledir, subdir)

if ~exist([dirname,subdir],'dir')
    mkdir([dirname,subdir]);
end

% get list of speech files
files = rdir([sampledir, '_*.mat']);
if length(files) == 0, error(['no file found in: ', sampledir]); end

% initialize storage variables
diglabels = {'one','two','three','four','five','six','seven','eight','nine','zero','oh'};
writename = cell(size(files));

% open mlf files
fid = fopen([dirname,subdir,'labels.mlf'],'w');
fprintf(fid,'#!MLF!#\n');
fid_sp = fopen([dirname,subdir,'labels_sp.mlf'],'w');
fprintf(fid_sp,'#!MLF!#\n');

% loop through files 
for samp = 1:length(files)
    
    writename{samp} = ['samp',num2str(samp)]; % This will be the file name
    
    % load data
    [~,~,~,data,~,annot]=load_data(files(samp).name,sampledir,strf.param.periph,...
        strf.param.normType,strf.param.int,strf.param.delay);
    data = data(:,strf.param.delay:end);
    disp(files(samp).name)
    
    % --- add stuff to mlf file
    fprintf(fid,['"*/',writename{samp},'.lab"\n']);
    fprintf(fid_sp,['"*/',writename{samp},'.lab"\n']);
    % put silence at the beginning and end of each transcription
    fprintf(fid,['sil','\n']);
    fprintf(fid_sp,['sil','\n']);
    % get rid of silence transcriptions
    transcriptions = annot.word2dig(annot.word2dig~=0);
    % put sp after all but the last transcription
    for i = 1:length(transcriptions)-1
        fprintf(fid,[diglabels{transcriptions(i)},'\n']);
        fprintf(fid_sp,[diglabels{transcriptions(i)},'\n']);
        fprintf(fid_sp,['sp','\n']);
    end
    % write the last transcription
    fprintf(fid,[diglabels{transcriptions(end)},'\n']);
    fprintf(fid_sp,[diglabels{transcriptions(end)},'\n']);
    fprintf(fid,['sil','\n']);
    fprintf(fid_sp,['sil','\n']);
    fprintf(fid_sp,'.\n');
    fprintf(fid,'.\n');
    
    % write the file.
    response = strf.H'*data;
    writesample(dirname,subdir,writename{samp},response,strf.param.int);
    
end

fclose(fid);
fclose(fid_sp);

% write data list
fid=fopen([dirname,subdir,'list.lst'],'w');
for i=1:length(files)
    fprintf(fid,[writename{i},'\n']);
end
fclose(fid);

end

function writesample(dirname,subdir,writename,response,int)

fid = fopen([dirname,subdir,writename],'wb');

% write number of samples (number of time points)
nsamp = size(response,2);
fwrite(fid, nsamp, 'integer*4',0,'b');

% write sample period (1)
period = int*10^4; % units of 100ns
fwrite(fid,period,'integer*4',0,'b');

% write sample size  = num channels*4byte floats
nbytes = size(response,1)*4;
fwrite(fid,nbytes,'integer*2',0,'b');

% write data type
type = 9;
fwrite(fid,type,'integer*2',0,'b');

% write data
fwrite(fid, response, 'real*4',0,'b');

fclose(fid);

end