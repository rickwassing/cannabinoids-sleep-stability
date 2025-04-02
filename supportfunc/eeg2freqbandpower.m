function EEG = eeg2freqbandpower(EEG, freqband)
% -------------------------------------------------------------------------
% Constants
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Apply wavelet to EEG data with a frequency step of:
freq_upper_limit = 0.5;
wavelet_freqstep = 0.1;
% Resample infraslow power modulation timeseries to a sampling rate of:
srate_resamp = 10; % Hz
% Low-pass filter settings
idx_bouts = round([1, EEG.event(strcmpi({EEG.event.type}, 'boundary')).latency, EEG.pnts+1]);
totit = EEG.nbchan * (length(idx_bouts)-1);
rt = now(); %#ok<TNOW1>
% -------------------------------------------------------------------------
% For each bout...
for i = 1:length(idx_bouts)-1
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Extract the bout
    this_idx = idx_bouts(i):(idx_bouts(i+1)-1);
    bout = EEG.data(:, this_idx)';
    for j = 1:EEG.nbchan
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        rt = remainingTime(rt, totit);
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        % Apply wavelet analysis on EEG data in a frequency band of interest
        pow_timeseries = morletgabortransform(bout(:, j), EEG.srate, freqband(1), freqband(2), wavelet_freqstep);
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        % average the power time course across all frequencies of interest
        EEG.data(j, this_idx) = mean(abs(pow_timeseries))';
    end
end


for j = 1:EEG.nbchan
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    bout = EEG.data(j, :)';
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    rt = remainingTime(rt, EEG.nbchan);
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Apply wavelet analysis on EEG data in a frequency band of interest
    pow_timeseries = morletgabortransform(bout, EEG.srate, freqband(1), freqband(2), wavelet_freqstep);
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % average the power time course across all frequencies of interest
    EEG.data(j, :) = mean(abs(pow_timeseries))';
end


% -------------------------------------------------------------------------
% Apply low-pass filter
EEG = pop_firws(EEG, ...
    'fcutoff', freq_upper_limit, ...
    'ftype', 'lowpass', ...
    'wtype', 'hamming', ...
    'forder', pop_firwsord('hamming', EEG.srate, 1), ...
    'minphase', 0);
% -------------------------------------------------------------------------
% De-mean
EEG.data = EEG.data - repmat(mean(EEG.data, 2), 1, EEG.pnts);
% -------------------------------------------------------------------------
% Resample
EEG = pop_resample(EEG, srate_resamp);
EEG.times = EEG.times ./1000;

end