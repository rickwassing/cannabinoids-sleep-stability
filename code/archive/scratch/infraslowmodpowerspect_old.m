function [POW, FREQ] = infraslowmodpowerspect_old(NREM, freqband)
% -------------------------------------------------------------------------
% Constants
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Number of FFT points to resample to
N = 2000;
% Upper limit of infraslow modulation frequency:
freq_upper_limit = 0.5;
% Apply wavelet to EEG data with a frequency step of:
wavelet_freqstep = 0.1;
% Resample infraslow power modulation timeseries to a sampling rate of:
srate_resamp = 10; % Hz
% Low-pass filter settings
filt_order = pop_firwsord('hamming', NREM.srate, 2*freq_upper_limit);
b_lowpass = fir1(filt_order, freq_upper_limit/(NREM.srate/2), 'low');
% -------------------------------------------------------------------------
% Extract bouts
idx_bouts = round([1, NREM.event(strcmpi({NREM.event.type}, 'boundary')).latency, NREM.pnts+1]);
bouts = {};
for i = 1:length(idx_bouts)-1
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Extract the bout
    this_idx = idx_bouts(i):(idx_bouts(i+1)-1);
    bouts{i} = NREM.data(21, this_idx)'; %#ok<AGROW>
end
% -------------------------------------------------------------------------
% Init output
FREQ = linspace(0, srate_resamp/2, N);
idx_f = find(FREQ <= freq_upper_limit);
POW = zeros(length(bouts), length(idx_f));
% -------------------------------------------------------------------------
% For each bout in the EEG data...
for i = 1:length(bouts)
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Apply wavelet analysis on EEG data in a frequency band of interest
    pow_timeseries = morletgabortransform(bouts{i}, NREM.srate, freqband(1), freqband(2), wavelet_freqstep);
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % average the power time course across all frequencies of interest
    pow_timeseries = mean(abs(pow_timeseries))';
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Apply low-pass filter
    pow_timeseries = filtfilt(b_lowpass, 1, double(pow_timeseries));
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Downsample and de-mean
    pow_timeseries = resample(pow_timeseries, srate_resamp, NREM.srate);
    pow_timeseries = pow_timeseries - mean(pow_timeseries);
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % The FFT frequencies of this bout
    fft_freqs = linspace(0, srate_resamp/2, (length(pow_timeseries)/2)+1);
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Apply FFT to the freq-band power timeseries
    pow_infraslow = (2*abs(fft(pow_timeseries))/length(pow_timeseries)).^2;
    pow_infraslow = pow_infraslow(1:length(fft_freqs));
    % Interpolate to standard length (to account for differences in bout length)
    pow_infraslow = interp1(fft_freqs, pow_infraslow, FREQ);
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Store in matrix
    POW(i,:) = pow_infraslow(idx_f);
end
% -------------------------------------------------------------------------
% Crop frequency range 
FREQ = FREQ(idx_f);
% -------------------------------------------------------------------------
% Calculate relative power (normalize to total power below 0.1 Hz)
POW = POW'; % transpose to frequencies by bouts
POW = POW ./ repmat(mean(abs(POW(FREQ <= 0.1, :))), size(POW, 1), 1);
POW = POW'; % transpose back to bouts by frequencies
POW = mean(POW, 1); %#ok<UDIM> % Average across bouts

end