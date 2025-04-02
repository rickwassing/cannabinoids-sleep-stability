function [m_Spindles, v_Duration, v_SpindleFreqs] = f_SpDetection_Ferrarelli(v_Signals, s_Fs, v_Hyp, varargin)
try
    if nargin < 4
        v_Aro = false(size(v_Signals));
    end
    if nargin < 3
        v_Hyp = -2 .* ones(size(v_Signals));
    end
    % Config
    if varargin{2} == 9999 % Pz
        b_DoPlot = true;
    else
        b_DoPlot = false;
    end
    s_LoRatio = 2;
    s_HiRatio = 4.5;
    v_HiPassFreq = [10, 30];
    v_BandPassFreq = [10, 16];
    v_Times = 1/s_Fs:1/s_Fs:length(v_Signals)/s_Fs;
    v_NremIdx = (v_Hyp == -2 | v_Hyp == -3); % & ~v_Aro;

    % Bandpass filter
    % Passband 11-15 Hz
    % Stopband 0-9 and 17-inf Hz
    forder = pop_firwsord('hamming', s_Fs, 2);
    bhi = fir1(forder, v_BandPassFreq/(s_Fs/2));
    v_BandPassFilt = filtfilt(bhi, 1, v_Signals);

    % Apply only a highpass to determine the spindle frequency
    bhi = fir1(forder, v_HiPassFreq/(s_Fs/2));
    v_HiPassFilt = filtfilt(bhi, 1, v_Signals);

    % Rectify signal
    v_Rectified = abs(v_BandPassFilt);

    % Create envelope from the peaks of rectified signal
    [v_Envelope, v_EnvelopeIdx] = findpeaks(v_Rectified);

    % Finds peaks of the envelope
    datader = diff(v_Envelope);
    posder = zeros(length(datader), 1);
    posder(datader > 0) = 1; % index of all points at which the rectified signal is increasing in amplitude
    diffder = diff(posder);
    v_EnvelopePeakIdx = v_EnvelopeIdx(find(diffder == -1) + 1); % peak index of Envelope signal
    v_EnvelopePeaks = v_Rectified(v_EnvelopePeakIdx); % peak amplitude of Envelope signal

    % Finds troughs of the envelope
    v_EnvelopeTroughIdx = v_EnvelopeIdx(find(diffder == 1) + 1); % trough index of Envelope signal
    v_EnvelopeTroughs = v_Rectified(v_EnvelopeTroughIdx); % peak trough of Envelope signal

    % Determine upper thresholds
    [v_Counts, v_Edges] = histcounts(v_EnvelopePeaks(v_NremIdx(v_EnvelopePeakIdx)), 120); % divide the distribution peaks of the Envelope signal in 120 bins
    [~, s_MaxCountIdx] = max(v_Counts); % select the most numerous bin
    s_MaxEnvAmp = v_Edges(s_MaxCountIdx) + mean(diff(v_Edges))/2; % peak of the amplitude distribution
    s_LoThres = s_LoRatio*s_MaxEnvAmp;

    % Determine lower threshold
    s_GoodSamplesThres = exp(mean(log(v_Rectified(v_NremIdx))) + 3*std(log(v_Rectified(v_NremIdx)))); % threshold this to deal with artefactual data samples
    s_HiThres = s_HiRatio * mean(v_Rectified(v_NremIdx & v_Rectified < s_GoodSamplesThres));

    % Define thresholds to find peaks in the envelope above the hi-threshold and where the troughs are below the lo-threshold
    v_TroughIdxBelowLoThres = v_EnvelopeTroughIdx(v_EnvelopeTroughs < s_LoThres);
    v_PeakIdxAboveHiThres = v_EnvelopePeakIdx(v_EnvelopePeaks > s_HiThres & v_NremIdx(v_EnvelopePeakIdx));

    % Initialize the spindle output
    m_Spindles = nan(2, length(v_PeakIdxAboveHiThres));
    v_SpindleFreqs = nan(1, length(v_PeakIdxAboveHiThres));
    % Initialize the while loop
    n = 0;
    i = 1;
    % Extract spindle onset and offsets, and their mean frequency
    while i <= length(v_PeakIdxAboveHiThres)
        % Get current peak, and the timepoints below the threshold (start and end samples)
        s_CurrentPeakIdx = v_PeakIdxAboveHiThres(i);
        s_TroughIdxBeforePeak = v_TroughIdxBelowLoThres(find(v_TroughIdxBelowLoThres > 1 & v_TroughIdxBelowLoThres < s_CurrentPeakIdx, 1, 'last'));
        s_TroughIdxAfterPeak  = v_TroughIdxBelowLoThres(find(v_TroughIdxBelowLoThres < length(v_Rectified) & v_TroughIdxBelowLoThres > s_CurrentPeakIdx, 1, 'first'));
        if isempty(s_TroughIdxBeforePeak) || isempty(s_TroughIdxAfterPeak)
            i = i+1;
            continue; % did not find a through on either side of the peak
        end
        % Store start and end samples in maxtrix
        m_Spindles(1, i) = s_TroughIdxBeforePeak;
        m_Spindles(2, i) = s_TroughIdxAfterPeak;
        % Calculate the average frequency
        [~, v_SpindlePeaks] = findpeaks(v_HiPassFilt(s_TroughIdxBeforePeak:s_TroughIdxAfterPeak));
        v_SpindleFreqs(i) = median(1./(diff(v_SpindlePeaks)./s_Fs));
        % Check how many peaks were in this interval, skip that number ahead
        s_NumPeaks = sum(v_PeakIdxAboveHiThres > s_TroughIdxBeforePeak & v_PeakIdxAboveHiThres < s_TroughIdxAfterPeak);
        i = i+s_NumPeaks;
    end

    % Check spindle duration
    v_Duration = m_Spindles(2, :) - m_Spindles(1, :);
    idx_rm = v_Duration < (0.3*s_Fs) | v_Duration > (3*s_Fs) | isnan(v_Duration);
    m_Spindles(:, idx_rm) = [];
    v_SpindleFreqs(idx_rm) = [];
    v_Duration(idx_rm) = [];

    if b_DoPlot
        % Apply broadband filter
        forder = pop_firwsord('hamming', s_Fs, 2);
        bhi = fir1(forder, [1, 30]/(s_Fs/2));
        v_BroadBandFilt = filtfilt(bhi, 1, v_Signals);
        % Initalize figure
        close all
        Fig = figure();
        Fig.Position = [1 200 1920 200];
        Fig.Color = [.61 .62 .64];
        Ax = axes(Fig);
        Ax.NextPlot = 'add';
        Ax.Color = Fig.Color;
        Ax.Layer = 'top';
        Ax.TickLength = [0, 0];
        Ax.Box = 'on';
        Ax.XGrid = 'on';
        Ax.GridColor = 'w';
        Ax.GridAlpha = 0.4;
        Ax.YTick = [];
        Ax.YLim = [-2, 30];
        Ax.XTick = unique(round(v_Times));
        Ax.XLim = [0, 30];
        % Plot
        plot(v_Times, v_Rectified, '-', 'Color', [0.9 .9 .9])
        plot(v_Times, (v_BroadBandFilt/10)+20, '-k')
        plot([0, max(v_Times)], [s_LoThres, s_LoThres], ':r')
        plot([0, max(v_Times)], [s_HiThres, s_HiThres], ':r')
        plot(v_Times(v_TroughIdxBelowLoThres), v_Rectified(v_TroughIdxBelowLoThres), '.k', 'MarkerSize', 3)
        plot(v_Times(v_PeakIdxAboveHiThres), v_Rectified(v_PeakIdxAboveHiThres), 'or', 'MarkerSize', 3)
        plot(v_Times(m_Spindles(:)), v_Rectified(m_Spindles(:)), 'or')
        text(v_Times(round(mean(m_Spindles))), ones(1, length(v_SpindleFreqs)), arrayfun(@(str) sprintf('%.1f', str), v_SpindleFreqs, 'UniformOutput', false), 'HorizontalAlignment', 'center')
        v_IsNrem = nan(size(v_NremIdx));
        v_IsNrem(v_NremIdx) = -1;
        plot(v_Times, v_IsNrem, '-', 'Color', [.3 .4 .5], 'LineWidth', 3);
        % Shift the x-axis limits to view each spindle and export to disk
        Ax.OuterPosition = [0.01, 0.01, 0.98, 0.98];
        for s = 1:size(m_Spindles, 2)
            Ax.XLim = [0, 30] + m_Spindles(1, s)/s_Fs - 15;
            kv = filename2struct(varargin{1}.setname);
            kv.desc = sprintf('spindle%i', s);
            kv.filetype = 'ferrarelli.png';
            filepath = sprintf('%s/images/spindledet/ferrarelli', varargin{1}.filepath);
            fullfilepath = sprintf('%s/%s', filepath, struct2filename(kv));
            if exist(filepath, 'dir') == 0
                mkdir(filepath)
            end
            exportgraphics(Fig, fullfilepath, 'Resolution', 144)
        end
    end
catch ME
    keyboard
end
end
