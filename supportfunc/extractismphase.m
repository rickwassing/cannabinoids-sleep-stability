function [EEG, NREM, SIGMA] = extractismphase(EEG, NREM, SIGMA, FIT)
% -------------------------------------------------------------------------
% Filter between ~0.015 and ~0.03 Hz, with a transition bandwith of 0.01
filtorder = pop_firwsord('hamming', SIGMA.srate, 0.01);
INFRA = SIGMA; % Copy structure to store filtered data
eeg_chans = find(strcmpi({SIGMA.chanlocs.type}, 'EEG')); % find EEG channels
% For each channel
rt = now; %#ok<TNOW1>
for i = 1:length(eeg_chans)
    % Extract the filter cutoff values and add padding of half the filter bandwidth
    fcutoff = [FIT(i).lower, FIT(i).upper] + [-0.005, 0.005];
    if fcutoff(1) < 0.0075
        fcutoff(1) = 0.0075;
    end
    % Select the channel and apply filter
    tmp = pop_firws(pop_select(SIGMA, 'channel', eeg_chans(i)), ...
        'fcutoff', fcutoff, ...
        'ftype', 'bandpass', ...
        'wtype', 'hamming', ...
        'forder', filtorder, ...
        'minphase', 0);
    INFRA.data(eeg_chans(i), :) = tmp.data - mean(tmp.data, 'omitnan'); % Store
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    rt = remainingTime(rt, length(eeg_chans), 'simple', true);
end
% -------------------------------------------------------------------------
% Replace NaN's with zero's
idx_nan = sum(isnan(INFRA.data), 1) > 0;
INFRA.data(:, idx_nan) = 0;
% -------------------------------------------------------------------------
% Apply Hilbert transform and get instantaneous amplitude and phase
AMPL = nan(length(eeg_chans), INFRA.pnts);
PHASE = nan(length(eeg_chans), INFRA.pnts);
for i = 1:length(eeg_chans)
    tmp = hilbert(INFRA.data(eeg_chans(i), :));
    AMPL(i, :) = abs(tmp)';
    PHASE(i, :) = angle(tmp)';
end
AMPL(:, idx_nan) = NaN;
PHASE(:, idx_nan) = NaN;
% -------------------------------------------------------------------------
% Store the amplitude and phase of each event
for i = 1:length(EEG.event)
    EEG.event(i).ism_phs = NaN;
    EEG.event(i).ism_amp = NaN;
    EEG.event(i).ism_pow = NaN;
    EEG.event(i).ism_phs_chan = NaN;
    EEG.event(i).ism_amp_chan = NaN;
    EEG.event(i).ism_pow_chan = NaN;
end
for i = 1:length(NREM.event)
    NREM.event(i).ism_phs = NaN;
    NREM.event(i).ism_amp = NaN;
    NREM.event(i).ism_pow = NaN;
    NREM.event(i).ism_phs_chan = NaN;
    NREM.event(i).ism_amp_chan = NaN;
    NREM.event(i).ism_pow_chan = NaN;
end
for i = 1:length(SIGMA.event)
    SIGMA.event(i).ism_phs = NaN;
    SIGMA.event(i).ism_amp = NaN;
    SIGMA.event(i).ism_pow = NaN;
    SIGMA.event(i).ism_phs_chan = NaN;
    SIGMA.event(i).ism_amp_chan = NaN;
    SIGMA.event(i).ism_pow_chan = NaN;
end
for i = 1:length(SIGMA.event)
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    idx = round(SIGMA.event(i).latency);
    SIGMA.event(i).ism_phs = circ_mean(PHASE(:, idx));
    SIGMA.event(i).ism_pow = mean(AMPL(:, idx), 'omitnan');
    SIGMA.event(i).ism_amp = mean(INFRA.data(:, idx), 'omitnan');
    SIGMA.event(i).ism_phs_chan = {PHASE(:, idx)};
    SIGMA.event(i).ism_pow_chan = {AMPL(:, idx)};
    SIGMA.event(i).ism_amp_chan = {INFRA.data(:, idx)};
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Add the same info to the NREM and EEG dataset
    idx_ev = find([EEG.event.id] == SIGMA.event(i).id);
    if ~isempty(idx_ev)
        EEG.event(idx_ev).ism_phs = SIGMA.event(i).ism_phs;
        EEG.event(idx_ev).ism_amp = SIGMA.event(i).ism_amp;
        EEG.event(idx_ev).ism_pow = SIGMA.event(i).ism_pow;
        EEG.event(idx_ev).ism_phs_chan = SIGMA.event(i).ism_phs_chan;
        EEG.event(idx_ev).ism_amp_chan = SIGMA.event(i).ism_amp_chan;
        EEG.event(idx_ev).ism_pow_chan = SIGMA.event(i).ism_pow_chan;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    idx_ev = find([NREM.event.id] == SIGMA.event(i).id);
    if ~isempty(idx_ev)
        NREM.event(idx_ev).ism_phs = SIGMA.event(i).ism_phs;
        NREM.event(idx_ev).ism_amp = SIGMA.event(i).ism_amp;
        NREM.event(idx_ev).ism_pow = SIGMA.event(i).ism_pow;
        NREM.event(idx_ev).ism_phs_chan = SIGMA.event(i).ism_phs_chan;
        NREM.event(idx_ev).ism_amp_chan = SIGMA.event(i).ism_amp_chan;
        NREM.event(idx_ev).ism_pow_chan = SIGMA.event(i).ism_pow_chan;
    end
end
end