function css_getsigmapowerusingwavelet(filepath, cfg)
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
% Get sigma power fluctuations: wavelet analysis in the freq-range of 11 to 15 Hz.
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Calculate Sigma power fluctuations (takes long time)
Settings = waveletsettings(EEG);
SIGMA = ft_wavelettransform(EEG, Settings);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Save
[SIGMA.filepath, SIGMA.setname] = fileparts(cfg.outfilepath);
SIGMA.filename = [SIGMA.setname, '.set'];
SaveDataset(SIGMA, 'all');

end
