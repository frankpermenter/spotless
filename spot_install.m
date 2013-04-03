% installation script for SPOT

potdir=pwd;
n=length(potdir);
s=filesep;                    % slash character
if ~strcmp('spot',potdir(n-3:n))||((s~='\')&&(s~='/')), 
    %    error('Please install SPOT in a "spot" directory!')
end
fprintf('\n Installing SPOT in %s:\n updating the path...',potdir)
addpath(potdir);
addpath([potdir s 'bin']);
addpath([potdir s 'util']);
addpath([potdir s 'mint']);
addpath([potdir s 'mss']);
addpath([potdir s 'spotopt']);
fprintf('\n compiling the binaries...')
cd('bin');
mex mss_gset.c 
mex mss_gsum.c
cd('..');
fprintf('\n Done.\n')
