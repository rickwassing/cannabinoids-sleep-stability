function css_getspecpowerusingwavelet(filepath, cfg)
% -------------------------------------------------------------------------
% Check if the 'filepath' variable is a string (path to file) or a stuct
% (preloaded data in EEGLAB structure)
if isstruct(filepath)
    EEG = filepath; % 'EEG' struct was used instead of path
    clear filepath;
else
    EEG = LoadDataset(filepath, 'all');
end
% -------------------------------------------------------------------------
% Get spectral power fluctuations: wavelet analysis in the specified freq-range.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Calculate Sigma power fluctuations (takes long time)
Settings = waveletsettings(EEG);
Settings.SpectrogramSettings.MinFreq = cfg.MinFreq;
Settings.SpectrogramSettings.MaxFreq = cfg.MaxFreq+0.01;
Settings.SpectrogramSettings.Cycles = cfg.Cycles;
POW = ft_wavelettransform(EEG, Settings);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Save
[POW.filepath, POW.setname] = fileparts(cfg.outfilepath);
POW.filename = [POW.setname, '.set'];
SaveDataset(POW, 'all');

end
