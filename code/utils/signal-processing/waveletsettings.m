function Settings = waveletsettings(EEG)
Settings.SpectrogramSettings.ChanSel = {EEG.chanlocs(strcmpi({EEG.chanlocs.type}, 'EEG')).labels};
Settings.SpectrogramSettings.SpectrogramType = 'mtmconvol';
Settings.SpectrogramSettings.TimeStep = 0.1;
Settings.SpectrogramSettings.FreqStep = 0.1; % 0.1
Settings.SpectrogramSettings.MinFreq = 11; % 11
Settings.SpectrogramSettings.MaxFreq = 15.01; % 15.01
Settings.SpectrogramSettings.Cycles = 16;
Settings.SpectrogramSettings.HasBaseline = false;
Settings.SpectrogramSettings.DoBaselineNorm = false;
Settings.SpectrogramSettings.BaselineMethod = 'absolute';
Settings.SpectrogramSettings.BaselineInterval = [-1, 0];
Settings.DoSpectrogram = true;
end