% MATLAB script to design a bandpass filter using a Butterworth IIR filter
% and plot the attenuation (magnitude response in dB) for different filter orders


% Define sampling frequency and cutoff frequencies
fs = 10; % Sampling frequency in Hz
f1 = 0.0138; % Lower cutoff frequency in Hz
f2 = 0.0259; % Upper cutoff frequency in Hz

% Define filter orders to test
filter_orders = [1, 2, 3, 4];

% Frequency vector for plotting
figure; hold on;

for N = filter_orders
    % Design a bandpass Butterworth filter
    [b, a] = butter(N, [f1 f2] / (fs/2), 'bandpass');
    
    % Compute and plot the magnitude response in dB
    [H, W] = freqz(b, a, 1024, fs);
    plot(W, 20*log10(abs(H)), 'DisplayName', ['Order ', num2str(N)]);
end

% Configure plot
xlabel('Frequency (Hz)');
ylabel('Attenuation (dB)');
title('Magnitude Response (dB) of Butterworth Bandpass Filter');
legend; grid on;

%[b, a] = butter(order, [f1 f2] / (fs/2), 'bandpass');

% Apply zero-phase filtering to EEG data
%filtered_data = filtfilt(b, a, eeg_data);