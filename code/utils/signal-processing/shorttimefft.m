function [SFT, FOOOF, f, t] = shorttimefft(data, times, fs)

% Define parameters
nbchan = size(data, 1);
ntrials = size(data, 3);
nfft = round(1.5 .* fs);
noverlap = round(1 .* fs);
window = hann(nfft);

rt = now(); %#ok<TNOW1> % Keep track of how much time is remaining
FOOOF = [];
for i = 1:ntrials
    rt = remainingTime(rt, ntrials);
    % Extact the data for this trial and transpose to time x chans
    tdata = squeeze(data(:, :, i))';
    % Apply short-time fourier transform
    [s, f, t] = stft(tdata, fs, ...
        'Window', window,...
        'OverlapLength', noverlap, ...
        'FFTLength', nfft, ...
        'FrequencyRange', 'onesided');
    % Reorder so channels are first dim, then freqs, then time
    s = permute(s, [3, 1, 2]);
    % Crop to freqs of interest
    s = s(:, f <= 25, :);
    f = double(f(f <= 25));
    % Set time vector
    t = t + times(1);
    % Convert to absolute
    s = abs(s);
    % Init output
    if i == 1
        SFT = nan(nbchan, length(f), length(t), ntrials, 'single');
    end
    % Store data in output
    SFT(:, :, :, i) = s;
    % Remove 1/f of spectrum averaged across time prior to arousal
    s2foof = squeeze(mean(s(:, :, t < 0), 3));
    FOOOF(i).f = applyfooof(s2foof, f); %#ok<AGROW>
end
end
