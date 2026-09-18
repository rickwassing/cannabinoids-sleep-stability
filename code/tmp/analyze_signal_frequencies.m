function [freq_analysis] = analyze_signal_frequencies(data, sampling_rate, plot_results)
% Analyze frequency content of EEG data to guide smoothing parameters
%
% Inputs:
%   data: [channels x samples] matrix (e.g., 18 x 600)
%   sampling_rate: Hz (e.g., if 600 samples span 60s, then fs = 10 Hz)
%   plot_results: boolean, whether to plot results (default: true)
%
% Outputs:
%   freq_analysis: struct with frequency analysis results

if nargin < 3
    plot_results = true;
end

[n_channels, n_samples] = size(data);
freq_analysis = struct();

% Calculate frequency vector
freq_vector = (0:n_samples-1) * sampling_rate / n_samples;
nyquist_freq = sampling_rate / 2;

% Initialize storage
all_power_spectra = zeros(n_channels, n_samples);
mean_power_spectrum = zeros(1, n_samples);

% Calculate power spectral density for each channel
for ch = 1:n_channels
    
    % Remove DC component and detrend
    signal = detrend(data(ch, :));
    
    % Apply window to reduce spectral leakage
    windowed_signal = signal .* hanning(n_samples)';
    
    % Calculate FFT
    fft_signal = fft(windowed_signal);
    
    % Calculate power spectral density
    power_spectrum = abs(fft_signal).^2 / (sampling_rate * n_samples);
    
    % For real signals, double the power (except DC and Nyquist)
    power_spectrum(2:end-1) = 2 * power_spectrum(2:end-1);
    
    all_power_spectra(ch, :) = power_spectrum;
end

% Average across channels
mean_power_spectrum = mean(all_power_spectra, 1);

% Only keep positive frequencies up to Nyquist
freq_positive = freq_vector(1:floor(n_samples/2)+1);
power_positive = mean_power_spectrum(1:floor(n_samples/2)+1);

% Calculate key frequency metrics
[max_power, max_idx] = max(power_positive(2:end)); % Exclude DC
dominant_freq = freq_positive(max_idx + 1); % +1 because we excluded DC

% Find frequency ranges with significant power
power_db = 10*log10(power_positive);
power_threshold = max(power_db) - 20; % -20dB from peak
significant_freqs = freq_positive(power_db > power_threshold);

% EMPIRICAL NOISE FLOOR DETECTION
% Method 1: Find where power spectrum flattens (derivative approach)
power_gradient = diff(power_db);
smooth_gradient = smoothdata(power_gradient, 'movmean', 5);
noise_start_idx = find(abs(smooth_gradient) < 1, 1, 'first'); % Where slope becomes flat
if isempty(noise_start_idx)
    noise_start_freq = freq_positive(end);
else
    noise_start_freq = freq_positive(noise_start_idx);
end

% Method 2: Find frequency where power drops to noise floor level
% Assume noise floor is the mean power in top 25% frequencies
noise_floor_power = mean(power_positive(round(0.75*length(power_positive)):end));
signal_cutoff_idx = find(power_positive <= noise_floor_power * 2, 1, 'first');
if isempty(signal_cutoff_idx)
    signal_cutoff_freq = freq_positive(end);
else
    signal_cutoff_freq = freq_positive(signal_cutoff_idx);
end

% Method 3: Find frequency where cumulative power reaches 95% of total
cumulative_power = cumsum(power_positive) / sum(power_positive);
power_95_idx = find(cumulative_power >= 0.95, 1, 'first');
if isempty(power_95_idx)
    power_95_freq = freq_positive(end);
else
    power_95_freq = freq_positive(power_95_idx);
end

% Take conservative estimate (lowest of the three methods)
estimates = [noise_start_freq, signal_cutoff_freq, power_95_freq];
estimates = estimates(estimates > 0);
empirical_signal_cutoff = 0.2; %max(estimates)/2;

% Calculate spectral centroid (weighted mean frequency)
spectral_centroid = sum(freq_positive .* power_positive) / sum(power_positive);

% Store results
freq_analysis.frequencies = freq_positive;
freq_analysis.power_spectrum = power_positive;
freq_analysis.power_db = power_db;
freq_analysis.dominant_frequency = dominant_freq;
freq_analysis.spectral_centroid = spectral_centroid;
freq_analysis.significant_freq_range = [min(significant_freqs), max(significant_freqs)];
freq_analysis.sampling_rate = sampling_rate;
freq_analysis.nyquist_frequency = nyquist_freq;

% Calculate empirical smoothing parameters based on noise floor
% Only smooth frequencies above the empirical signal cutoff
if empirical_signal_cutoff < nyquist_freq
    empirical_smoothing_samples = max(1, round(sampling_rate / (2 * empirical_signal_cutoff)));
else
    empirical_smoothing_samples = 1; % No smoothing needed
end

freq_analysis.noise_start_frequency = noise_start_freq;
freq_analysis.signal_cutoff_frequency = signal_cutoff_freq;
freq_analysis.power_95_frequency = power_95_freq;
freq_analysis.empirical_signal_cutoff = empirical_signal_cutoff;
freq_analysis.empirical_smoothing_span = empirical_smoothing_samples;

if plot_results
    figure('Position', [100, 100, 1000, 600]);
    
    % Plot 1: Power spectrum
    subplot(2,3,1);
    semilogy(freq_positive, power_positive);
    hold on;
    semilogy(dominant_freq, max_power, 'ro', 'MarkerSize', 8, 'LineWidth', 2);
    xlabel('Frequency (Hz)');
    ylabel('Power');
    title('Average Power Spectrum');
    grid on;
    legend('Power Spectrum', sprintf('Dominant: %.3f Hz', dominant_freq), 'Location', 'best');
    
    % Plot 2: Power spectrum in dB
    tmpax = subplot(2,3,2);
    plot(freq_positive, power_db);
    hold on;
    plot([min(freq_positive), max(freq_positive)], [power_threshold, power_threshold], 'r--');
    xlabel('Frequency (Hz)');
    ylabel('Power (dB)');
    title('Power Spectrum (dB scale)');
    grid on;
    legend('Power (dB)', '-20dB threshold', 'Location', 'best');
    tmpax.XScale = 'log';
    
    % Plot 3: Individual channel spectra
    subplot(2,3,3);
    for ch = 1:min(n_channels, 6) % Plot max 6 channels for clarity
        power_ch = all_power_spectra(ch, 1:floor(n_samples/2)+1);
        semilogy(freq_positive, power_ch, 'Color', [0.7, 0.7, 0.7]);
        hold on;
    end
    semilogy(freq_positive, power_positive, 'k-', 'LineWidth', 2);
    xlabel('Frequency (Hz)');
    ylabel('Power');
    title('Individual Channels + Average');
    grid on;
    
    % Plot 4: Example smoothing comparison
    subplot(2,3,4:6);
    example_channel = 1;
    original_signal = data(example_channel, :);
    time_vector = (0:n_samples-1) / sampling_rate;
    
    plot(time_vector, original_signal, 'b-', 'LineWidth', 1);
    hold on;
    
    % Test different smoothing spans based on empirical analysis
    smoothing_spans = [1, empirical_smoothing_samples];
    colors = {'b', 'r', 'g'};
    labels = {'Original'};
    
    for i = 1:length(smoothing_spans)
        if smoothing_spans(i) > 1
            smoothed = smoothdata(original_signal, 'movmean', smoothing_spans(i));
            plot(time_vector, smoothed, colors{i}, 'LineWidth', 1.5);
            if i == 1
                labels{end+1} = sprintf('Empirical (span=%d)', smoothing_spans(i));
            else
                labels{end+1} = sprintf('Heavy (span=%d)', smoothing_spans(i));
            end
        end
    end
    
    xlabel('Time (s)');
    ylabel('Amplitude');
    title(sprintf('Smoothing Comparison (Channel %d)', example_channel));
    legend(labels, 'Location', 'best');
    grid on;
    
    % Display summary
    fprintf('\n=== EMPIRICAL FREQUENCY ANALYSIS ===\n');
    fprintf('Sampling Rate: %.2f Hz\n', sampling_rate);
    fprintf('Dominant Frequency: %.3f Hz\n', dominant_freq);
    fprintf('Spectral Centroid: %.3f Hz\n', spectral_centroid);
    fprintf('Significant Freq Range: %.3f - %.3f Hz\n', freq_analysis.significant_freq_range);
    fprintf('\n--- EMPIRICAL NOISE ANALYSIS ---\n');
    fprintf('Noise floor starts at: %.3f Hz (gradient method)\n', noise_start_freq);
    fprintf('Signal cutoff at: %.3f Hz (noise floor method)\n', signal_cutoff_freq);
    fprintf('95%% power contained below: %.3f Hz\n', power_95_freq);
    fprintf('Conservative signal cutoff: %.3f Hz\n', empirical_signal_cutoff);
    fprintf('Empirical smoothing span: %d samples\n', empirical_smoothing_samples);
    fprintf('===================================\n\n');
end

end

% Example usage:
%
% % If your 600 samples span 60 seconds:
% sampling_rate = 10; % Hz (600 samples / 60 seconds)
% freq_analysis = analyze_signal_frequencies(sigma_data, sampling_rate);
%
% % Use the recommended smoothing span:
% optimal_span = freq_analysis.recommended_smoothing_span;
% smoothed_data = smoothdata(sigma_data, 2, 'movmean', optimal_span);
%
% % Or test different smoothing levels:
% conservative_span = round(optimal_span * 0.5); % Less smoothing
% aggressive_span = round(optimal_span * 2);     % More smoothing

% Alternative: Welch's method for more robust PSD estimation
function [freq_analysis_welch] = analyze_frequencies_welch(data, sampling_rate)
% More robust frequency analysis using Welch's method

[n_channels, n_samples] = size(data);

% Parameters for Welch's method
window_length = min(256, floor(n_samples/4));
overlap = floor(window_length/2);

all_psds = [];

for ch = 1:n_channels
    [psd, freq] = pwelch(data(ch, :), window_length, overlap, [], sampling_rate);
    all_psds = [all_psds; psd'];
end

mean_psd = mean(all_psds, 1);
freq_analysis_welch.frequencies = freq;
freq_analysis_welch.power_spectrum = mean_psd;
freq_analysis_welch.method = 'Welch';

end