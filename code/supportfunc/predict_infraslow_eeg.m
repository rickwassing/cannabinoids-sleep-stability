function [predicted_signal, predicted_complex, analytic_signal] = predict_infraslow_eeg(signal, varargin)
% PREDICT_INFRASLOW_EEG Predicts future samples of bandpass filtered EEG sigma power
% using phase-amplitude decomposition
%
% Inputs:
%   signal - Input vector of filtered EEG data (should be ~60 seconds at 10 Hz)
%   varargin - Optional parameters:
%              'sampling_rate' - Sampling rate in Hz (default: 10)
%              'prediction_length' - Number of samples to predict (default: 50)
%              'do_plot' - Whether to plot results (default: false)
%
% Output:
%   predicted_signal - Predicted signal for the next prediction_length samples

% Parse input arguments
p = inputParser;
addRequired(p, 'signal', @isnumeric);
addParameter(p, 'sampling_rate', 10, @isnumeric);
addParameter(p, 'prediction_length', 50, @isnumeric);
addParameter(p, 'do_plot', false, @islogical);
addParameter(p, 'padding_samples', 0, @isnumeric)
parse(p, signal, varargin{:});

fs = p.Results.sampling_rate;
pred_len = p.Results.prediction_length;
do_plot = p.Results.do_plot;
padding_samples = p.Results.padding_samples;

% Handle padding for better Hilbert transform edge behavior
if padding_samples > 0
    % Signal comes in with padding still attached
    signal_full = signal(:);

    % Extract the actual signal portion (middle)
    signal_actual = signal_full(padding_samples+1:end-padding_samples);
    N = length(signal_actual);

    % Compute Hilbert transform on full padded signal for better edges
    analytic_signal_full = hilbert(signal_full);

    % Extract the middle portion for analysis (no edge artifacts)
    analytic_signal = analytic_signal_full(padding_samples+1:end-padding_samples);

    % Use actual signal for subsequent processing
    signal = signal_actual;
else
    % Original behavior - no padding
    N = length(signal);
    analytic_signal = hilbert(signal);
end

% Extract instantaneous phase and amplitude
inst_phase = unwrap(angle(analytic_signal));
inst_amplitude = abs(analytic_signal);

% === PHASE PREDICTION (Weighted Linear Regression) ===
% Fit weighted linear model to unwrapped phase
t = (0:N-1)' / fs;

% Create exponentially decaying weights (recent samples weighted more)
decay_factor = 0.98; % Adjust between 0.8-0.99 (higher = more emphasis on recent)
weights = decay_factor.^((N-1):-1:0)'; % More weight to recent samples
weights = weights / sum(weights); % Normalize

% Weighted linear regression for phase (y = a + b*t)
X = [ones(N,1), t]; % Design matrix for linear fit
W = diag(weights); % Weight matrix
phase_coeffs = (X'*W*X) \ (X'*W*inst_phase); % Weighted least squares

inst_freq = phase_coeffs(2) / (2*pi); % Convert slope to Hz

% Predict future phase linearly
t_future = (N:N+pred_len-1)' / fs;
predicted_phase = phase_coeffs(1) + phase_coeffs(2) * t_future;

% === AMPLITUDE PREDICTION (Nonlinear) ===
% Use a combination of methods for amplitude prediction

% Method 1: Exponential smoothing with trend
alpha = 0.2; % Smoothing parameter
beta = 0.1;  % Trend parameter

% Initialize
S = inst_amplitude(1);
T = inst_amplitude(2) - inst_amplitude(1);
smoothed_amp = zeros(N, 1);
trends = zeros(N, 1);

smoothed_amp(1) = S;
trends(1) = T;

% Apply Holt's exponential smoothing
for i = 2:N
    S_prev = S;
    T_prev = T;
    S = alpha * inst_amplitude(i) + (1 - alpha) * (S_prev + T_prev);
    T = beta * (S - S_prev) + (1 - beta) * T_prev;
    smoothed_amp(i) = S;
    trends(i) = T;
end

% Method 2: Autoregressive component for residuals
residuals = inst_amplitude - smoothed_amp;

% Fit AR model to residuals (order selection based on signal length)
ar_order = min(10, floor(N/10)); % Adaptive order selection
try
    ar_model = ar(residuals, ar_order, 'ls');
    ar_coeffs = ar_model.A(2:end); % Exclude the first coefficient (which is 1)
    ar_coeffs = ar_coeffs(:); % Ensure column vector
catch
    % Fallback to simpler approach if AR fitting fails
    ar_order = min(3, floor(N/20));
    ar_coeffs = lpc(residuals, ar_order);
    ar_coeffs = ar_coeffs(2:end);
    ar_coeffs = ar_coeffs(:); % Ensure column vector
end

% Predict amplitude
predicted_amplitude = zeros(pred_len, 1);

% Extend trend component
final_trend = trends(end);
final_level = smoothed_amp(end);

% Initialize residual prediction buffer (ensure column vector)
residual_buffer = residuals(end-ar_order+1:end);
residual_buffer = residual_buffer(:); % Ensure column vector

for i = 1:pred_len
    % Trend-based prediction
    trend_pred = final_level + i * final_trend;
    
    % AR prediction with damping
    if length(ar_coeffs) > 0
        residual_pred = -sum(ar_coeffs .* residual_buffer);
        % Add exponential damping to AR component
        damping_factor = 0.95^i; % Decay AR influence over time
        residual_pred = residual_pred * damping_factor;
        
        % Update buffer
        residual_buffer = [residual_buffer(2:end); residual_pred];
    else
        residual_pred = 0;
    end
    
    predicted_amplitude(i) = max(0, trend_pred + residual_pred);
end

% === RECONSTRUCT SIGNAL ===
% Convert back to complex form and take real part
predicted_complex = predicted_amplitude .* exp(1j * predicted_phase);
predicted_signal = real(predicted_complex);

% === PLOTTING ===
if do_plot
    Fig = figure;
    Fig.Position = [250 250 1050 700];

    % Time vectors
    t_orig = (0:N-1) / fs;
    t_pred = (N:N+pred_len-1) / fs;
    t_all = [t_orig, t_pred];
    signal_all = [asrow(signal), asrow(predicted_signal)];

    clear Ax
    Ax(1) = subplot(3,1,1);
    plot(t_orig, signal, 'b-', 'LineWidth', 1.5); hold on;
    plot(t_pred, predicted_signal, 'r-', 'LineWidth', 1.5);
    plot([t_orig(end), t_pred(1)], [signal(end), predicted_signal(1)], 'g:', 'LineWidth', 2);
    xlabel('Time (s)');
    ylabel('Amplitude');
    title('Original vs Predicted Signal');

    grid on;

    Ax(2) = subplot(3,1,2);
    plot(t_orig, mod(inst_phase, 2*pi), 'b-', 'LineWidth', 1); hold on;
    plot(t_pred, mod(predicted_phase, 2*pi), 'r-', 'LineWidth', 1);
    xlabel('Time (s)');
    ylabel('Phase (rad)');
    title(sprintf('Instantaneous Phase (Est. Freq: %.4f Hz)', inst_freq));
    Ax(2).YLim = [-0.1, 6.25];
    grid on;

    Ax(3) = subplot(3,1,3);
    plot(t_orig, inst_amplitude, 'b-', 'LineWidth', 1); hold on;
    plot(t_pred, predicted_amplitude, 'r-', 'LineWidth', 1);
    plot(t_orig, smoothed_amp, 'g:', 'LineWidth', 1);
    xlabel('Time (s)');
    ylabel('Amplitude');
    title('Instantaneous Amplitude');

    grid on;

    linkaxes(Ax, 'x');

    Ax(1).XLim = [t_all(end)-50, t_all(end)+5];
    sgtitle('EEG Infraslow Signal Prediction via Phase-Amplitude Decomposition');
end

end