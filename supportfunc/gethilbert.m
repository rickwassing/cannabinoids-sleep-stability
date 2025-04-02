function [POW, bouts] = gethilbert(SPEC, HYP, freqs)

% Select NREM stage 2 and 3 time segments during sws and asc sleep
bouts = getbouts(HYP, SPEC, 120);

POW = {};
if ischar(freqs)
    for i = 1:size(bouts, 1)
        idx_time = SPEC.spectimes >= bouts(i, 1) & SPEC.spectimes <= bouts(i, 2);
        POW{i} = SPEC.smooth.(freqs)(idx_time); %#ok<AGROW> 
    end
else
    idx_freq = SPEC.specfreqs >= freqs(1) & SPEC.specfreqs <= freqs(2);
    for i = 1:size(bouts, 1)
        idx_time = SPEC.spectimes >= bouts(i, 1) & SPEC.spectimes <= bouts(i, 2);
        POW{i} = squeeze(mean(SPEC.specdata(:, idx_freq, idx_time), 2, 'omitnan')); %#ok<AGROW> 
    end
end

% rmidx = diff(point') < 300*EEG.srate;
% point(rmidx, :) = [];
% channel = find(ismember({EEG.chanlocs.labels}, [ChansFz, ChansCz, ChansPz]));
% EEG = pop_select(EEG, 'channel', channel, 'point', point);
% EEG.spectimes = [];
% EEG.specfreqs = [];
% EEG.specchans = [];
% EEG.specdata = [];
% EEG.times = EEG.xmin:1/EEG.srate:EEG.xmax;
% 
% % -------------------------------------------------------------------------
% % Filter in sigma freq range (spindle range = 11 - 15 Hz; Purcell 2017)
% Settings.DoFilter = true;
% Settings.FilterSettings.DoNotch = false;
% Settings.FilterSettings.DoBandpass = true;
% Settings.FilterSettings.Highpass = freqs(1);
% Settings.FilterSettings.Lowpass = freqs(2);
% Settings.FilterSettings.WindowType = 'hamming';
% Settings.FilterSettings.TransitionBW = 1;
% EEG = Proc_TemporalFilter(EEG, Settings);
% % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% % Apply Hilbert
% for i = 1:EEG.nbchan
%     EEG.data(i, :) = abs(hilbert(EEG.data(i,:)));
% end
% 
% % -------------------------------------------------------------------------
% % Average topological locations (Fz, Cz, Pz)
% tmp = nan(3, EEG.pnts, 'single');
% tmp(1, :) = mean(EEG.data(ismember({EEG.chanlocs.labels}, ChansFz), :));
% tmp(2, :) = mean(EEG.data(ismember({EEG.chanlocs.labels}, ChansCz), :));
% tmp(3, :) = mean(EEG.data(ismember({EEG.chanlocs.labels}, ChansPz), :));
% EEG.data = tmp;
% EEG.nbchan = 3;
% EEG.chanlocs = EEG.chanlocs(ismember({EEG.chanlocs.labels}, {'E21', 'Cz', 'E101'}));
% EEG.chanlocs = [EEG.chanlocs(1), EEG.chanlocs(3), EEG.chanlocs(2)];

end