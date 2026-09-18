function EEG = css_preproc(EEG, cfg)
% -------------------------------------------------------------------------
% Load data if not already done
if ~isstruct(EEG)
    % assume its the filepath
    EEG = LoadDataset(EEG, 'all');
end
EEG.subject = '';
EEG.comments = '';
EEG.urchanlocs = [];
EEG.urevent = [];
EEG.specdata = [];
EEG.specfreqs = [];
EEG.spectimes = [];
EEG.specchans = [];
% -------------------------------------------------------------------------
% Keep only some events
allowedevents = {'wake', 'n1', 'n2', 'n3', 'rem', 'alpha', 'alphaemg', 'arousal', 'arousalemg', 'loff', 'lon', 'reject'};
idx_keep = ismember({EEG.event.type}, allowedevents);
EEG.event = EEG.event(idx_keep);
% -------------------------------------------------------------------------
% Extract ECG channel
idx_ecg = find(strcmpi({EEG.chanlocs.type}, 'ecg'));
ECG = pop_select(EEG, 'channel', idx_ecg);
HR = calcinsthr(ECG);
HR = pop_select(HR, 'channel', 2);
% -------------------------------------------------------------------------
% Extract EEG channels
cranial_chans = hdeeg_scalpchannels('egi257');
idx_eeg = find(ismember({EEG.chanlocs.labels}, cranial_chans));
EEG = pop_select(EEG, 'channel', idx_eeg);
% -------------------------------------------------------------------------
% Average reference and interpolating rejected channels
Settings = struct();
Settings.DoInterpolate = true;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if isfield(EEG.etc, 'rej_channels')
    rejchans = EEG.etc.rej_channels;
else
    rejchans = [];
end
EEG = Proc_AverageRef(EEG, [], rejchans);
EEG = Proc_RejectChannels(EEG, Settings, rejchans, []);
% -------------------------------------------------------------------------
% Apply broadband filter
fcutoff = [0.3, 48];
forder = pop_firwsord('hamming', EEG.srate, 0.4);
EEG = pop_firws(EEG, ...
    'fcutoff', fcutoff, ...
    'ftype', 'bandpass', ...
    'wtype', 'hamming', ...
    'forder', forder, ...
    'minphase', 0);
% -------------------------------------------------------------------------
% Resample
EEG = pop_resample(EEG, 128);
HR = pop_resample(HR, 128);
% -------------------------------------------------------------------------
% Make sure time is in seconds
EEG.times = EEG.xmin:1/EEG.srate:EEG.xmax;
ECG.times = ECG.xmin:1/ECG.srate:ECG.xmax;
HR.times = HR.xmin:1/HR.srate:HR.xmax;
% -------------------------------------------------------------------------
% Sidecar files
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
EEG.etc.JSON.TaskName = 'psg';
EEG.etc.JSON.EEGReference = 'averef';
EEG.etc.JSON.EEGChannelCount = 178;
EEG.etc.JSON.PNSChannelCount = 0;
EEG.etc.JSON.RecordingDuration = EEG.xmax;
EEG.etc.JSON.RecordingType = 'continuous';
EEG.etc.JSON.TrialCount = 1;
try
    EEG.etc.JSON.SoftwareFilters.BandPass.PassBand = [0.3, 48];
    EEG.etc.JSON.SoftwareFilters.BandPass.TransitionBandwidth = 0.4;
    EEG.etc.JSON.SoftwareFilters.BandPass.Order = forder;
catch
    % Its ok
end
EEG.etc.JSON.IsDownsampled = true;
EEG.etc.JSON.SamplingFrequency = 128;
EEG.etc.JSON.ECGChannelCount = 0;
EEG.etc.JSON.EMGChannelCount = 0;
EEG.etc.JSON.EOGChannelCount = 0;
EEG.etc.JSON.MiscChannelCount = 0;
if isfield(EEG.etc.JSON, 'SpectralAnalysis')
    EEG.etc.JSON = rmfield(EEG.etc.JSON, 'SpectralAnalysis');
end
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ECG.etc.JSON.TaskName = 'psg';
ECG.etc.JSON.EEGReference = 'averef';
ECG.etc.JSON.EEGChannelCount = 0;
ECG.etc.JSON.PNSChannelCount = 1;
ECG.etc.JSON.RecordingDuration = ECG.xmax;
ECG.etc.JSON.RecordingType = 'continuous';
ECG.etc.JSON.TrialCount = 1;
ECG.etc.JSON.ECGChannelCount = 1;
ECG.etc.JSON.EMGChannelCount = 0;
ECG.etc.JSON.EOGChannelCount = 0;
ECG.etc.JSON.MiscChannelCount = 0;
if isfield(ECG.etc.JSON, 'SpectralAnalysis')
    ECG.etc.JSON = rmfield(ECG.etc.JSON, 'SpectralAnalysis');
end
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
HR.etc.JSON.TaskName = 'psg';
HR.etc.JSON.EEGReference = 'averef';
HR.etc.JSON.EEGChannelCount = 0;
HR.etc.JSON.PNSChannelCount = 1;
HR.etc.JSON.RecordingDuration = HR.xmax;
HR.etc.JSON.RecordingType = 'continuous';
HR.etc.JSON.TrialCount = 1;
HR.etc.JSON.IsDownsampled = true;
HR.etc.JSON.SamplingFrequency = 128;
HR.etc.JSON.ECGChannelCount = 0;
HR.etc.JSON.EMGChannelCount = 0;
HR.etc.JSON.EOGChannelCount = 0;
HR.etc.JSON.MiscChannelCount = 1;
if isfield(HR.etc.JSON, 'SpectralAnalysis')
    HR.etc.JSON = rmfield(HR.etc.JSON, 'SpectralAnalysis');
end
% -------------------------------------------------------------------------
% Save
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% EEG
[EEG.filepath, EEG.setname] = fileparts(cfg.outfilepath);
EEG.filename = [EEG.setname, '.set'];
SaveDataset(EEG, 'all');
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% ECG
kv = filename2struct(EEG.setname);
kv.desc = [kv.desc, 'ecg'];
kv.filetype = 'ecg';
ECG.filepath = EEG.filepath;
ECG.setname = struct2filename(kv);
ECG.filename = [ECG.setname, '.set'];
SaveDataset(ECG, 'all');
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% HR
kv = filename2struct(EEG.setname);
kv.desc = [kv.desc, 'hr'];
kv.filetype = 'hr';
HR.filepath = EEG.filepath;
HR.setname = struct2filename(kv);
HR.filename = [HR.setname, '.set'];
SaveDataset(HR, 'all');

end