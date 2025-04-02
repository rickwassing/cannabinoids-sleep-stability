% This code reproduces the findings reported in 
% DOI: xxx.x.x.x

% Author: Rick Wassing (rick.wassing@woolcock.org.au)

% =========================================================================
% INITIALIZE
% -------------------------------------------------------------------------
% Initialize Matlab, add code to path
cd('S:\Sleep\3. ACTIVE STUDIES\CUPID\Arousal paper backup_CANSLEEP')
addpath(genpath('code-paper'))
css_init('superpc')

%%

% Specify which computers you have available
hosts = {...
    'Ricks-MacBook-Air.local', ...
    'mnc.local', ...
    'WIMR-SUPERMICRO', ...
    'wimr-HP-Z6'};
[~, this_host] = system('hostname');
idx_host = find(strcmpi(strip(this_host), hosts));


%% =========================================================================
% PROCESSING
% -------------------------------------------------------------------------
% Create a list of files to process
clc
Files = dir('derivatives/EEG-preproc/sub-*/ses-*/sub-*_desc-chaninterpavref_eeg.set');
errors = {};
% -------------------------------------------------------------------------
% PART 1: Extract 300 seconds of NREM bouts,
for i = 1:length(Files)
    try
        css_proc(fullfile(Files(i).folder, Files(i).name)); % TODO rename this function
    catch ME
        disp('****************************************************')
        disp(getReport(ME))
        disp('****************************************************')
        err.it = i;
        err.File = Files(i);
        err.ME = ME;
        errors = [errors; {err}]; %#ok<AGROW>
    end
end 

%% apply Wavelet analysis and
% estimate the frequency and amplitude of sigma power modulation.

%% -------------------------------------------------------------------------
% PART 2: For each arousal, get the pre-arousal spectrogram and distinguish
% between continued sleep versus awakening arousals.

% -------------------------------------------------------------------------
% PART 3: For each arousal in NREM-stage 2 sleep, extract a minimum of 100
% seconds of NREM sleep, extract sigma power and calculate instantaneous
% phase at the onset of the arousal. Distinguish between continued sleep
% versus awakening arousals.
% -------------------------------------------------------------------------
% PART 4: For each transition to REM sleep check the infraslow modulation
% of sigma power, and distinguish between short and long REM sleep bouts.

% =========================================================================
% ANALYSIS
% -------------------------------------------------------------------------
% Load individual hypnogram data
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Describe PSG sleep parameters
% -------------------------------------------------------------------------
% Load individual arousal data
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Describe arousal types (continued sleep vs state-shifts), when they occur, 




