function [FIT] = infraslowmodpowerspect(SIGMA, doPlot, doFilter)
% -------------------------------------------------------------------------
% Constants
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Number of FFT points to resample to
N = 50000;
% Upper limit of infraslow modulation frequency:
freq_upper_limit = 0.5;
% Make figure?
if nargin < 3
    doFilter = true;
end
if nargin < 2
    doPlot = true;
end
% What method to use to extract ISM peaks?
fitMethod = 'gaussianfit'; % Options 'gaussianfit' or 'smoothdata'
ngaussians = 2; % 1 or 2
% -------------------------------------------------------------------------
% Extract bouts
idx_bouts = round([1, SIGMA.event(strcmpi({SIGMA.event.type}, 'boundary')).latency, SIGMA.pnts+1]);
bouts = {};
for i = 1:length(idx_bouts)-1
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Extract the bout and demean the data
    this_idx = idx_bouts(i):(idx_bouts(i+1)-1);
    x = SIGMA.data(:, this_idx)';
    % Remove nan's
    x(sum(isnan(x), 2) > 0, :) = [];
    % High-pass
    if doFilter
        forder = pop_firwsord('hamming', SIGMA.srate, 0.01);
        if length(x) <= forder*3
            x = detrend(x);
        else
            b = fir1(forder, 0.005/(SIGMA.srate/2), 'high');
            x = filtfilt(b, 1, double(x));
        end
    else
        x = detrend(x, 0); % demean
    end
    bouts{i} = x; %#ok<AGROW>
end
% -------------------------------------------------------------------------
% Init output
nfft = 300.*SIGMA.srate; % 2*N-1
window = hanning(nfft);
overlap = round(nfft/2);
FREQ_FFT = linspace(0, SIGMA.srate/2, N);
idx_f_fft = find(FREQ_FFT <= freq_upper_limit);
FREQ_PWL = linspace(0, SIGMA.srate/2, nfft/2+1);
idx_f_pwl = find(FREQ_PWL <= freq_upper_limit);
% -------------------------------------------------------------------------
% First and last bouts may have NaN values because of the edge artefacts of
% the wavelet analysis (half the wavelet length). Remove bouts that are too
% short.
idx_rm = cellfun(@(b) length(b), bouts) < nfft;
bouts(idx_rm) = [];
POW = zeros(length(bouts), length(idx_f_fft), SIGMA.nbchan);
PXX = zeros(length(bouts), length(idx_f_pwl), SIGMA.nbchan);
% -------------------------------------------------------------------------
% For each bout in the EEG data...
for i = 1:length(bouts)
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % The FFT frequencies of this bout
    fft_freqs = linspace(0, SIGMA.srate/2, (size(bouts{i}, 1)/2)+1);
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Apply FFT to the freq-band power timeseries
    pow_infraslow = (2*abs(fft(bouts{i}))./size(bouts{i}, 1)).^2;
    pow_infraslow = pow_infraslow(1:length(fft_freqs), :);
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Apply PWELCH to the freq-band power timeseries
    [pwel_infraslow, pwel_freqs] = pwelch(bouts{i}, window, overlap, nfft, SIGMA.srate, 'power');
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Interpolate FFT to standard length (to account for differences in bout length)
    pow_infraslow = interp1(fft_freqs, pow_infraslow, FREQ_FFT);
    % make sure its a column vector
    if length(size(pow_infraslow)) == 2 && any(size(pow_infraslow) == 1)
        pow_infraslow = ascolumn(pow_infraslow);
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Store in matrix
    POW(i, :, :) = pow_infraslow(idx_f_fft, :, :);
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Store in matrix
    PXX(i, :, :) = pwel_infraslow(idx_f_pwl, :, :);
end
% -------------------------------------------------------------------------
% Crop frequency range
FREQ_FFT = FREQ_FFT(idx_f_fft);
FREQ_PWL = FREQ_PWL(idx_f_pwl);
% -------------------------------------------------------------------------
% Calculate relative power (normalize to total power below 0.1 Hz)
POW = squeeze(mean(POW, 1)); % Average across bouts
mu = mean(mean(abs(POW(FREQ_FFT > 0.0075 & FREQ_FFT <= 0.1, :)), 1));
POW = POW ./ mu;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
PXX = squeeze(mean(PXX, 1)); % Average across bouts
mu = mean(mean(abs(PXX(FREQ_PWL > 0.0075 & FREQ_PWL <= 0.1, :)), 1));
PXX = PXX ./ mu;
% -------------------------------------------------------------------------
% Select which measure to use (FFT or PWELCH)
MEAS = PXX;
FREQ = FREQ_PWL;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% We're fitting gaussians to frequencies under 0.1 Hz
idx_freqfit = FREQ <= 0.1;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Init the output structure
FIT = struct();
% Init figure
if doPlot
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    close all
    Fig = figure('Color', 'w', 'Position', [560, 640, 290 200]);
    Ax = axes('NextPlot', 'add', 'FontSize', 10, 'XTick', 0:0.01:0.1, 'YTick', 1:10, 'XGrid', 'on', 'TickLength', [0 0], 'Box', 'on', 'Color', [0.94 0.95 0.97]);
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    H.data = plot(Ax, NaN, NaN, '-', 'Color', [0.6 0.6 0.6]);
    H.fit = plot(Ax, NaN, NaN, '-r', 'LineWidth', 1);
    H.gaus(1) = plot(Ax, NaN, NaN, ':r', 'LineWidth', 0.5);
    H.gaus(2) = plot(Ax, NaN, NaN, ':r', 'LineWidth', 0.5);
    H.gaus(3) = plot(Ax, NaN, NaN, ':r', 'LineWidth', 0.5);
    H.gaus(4) = plot(Ax, NaN, NaN, ':r', 'LineWidth', 0.5);
    H.thres = plot(Ax, NaN, NaN, '--k', 'LineWidth', 0.75);
    H.lolim = plot(Ax, NaN, NaN, '--r', 'LineWidth', 0.75);
    H.hilim = plot(Ax, NaN, NaN, '--r', 'LineWidth', 0.75);
    H.peak = plot(Ax, NaN, NaN, 'or', 'MarkerFaceColor', 'r', 'MarkerSize', 6);
    H.span = plot(Ax, [NaN, NaN], [2, 2], '-k', 'LineWidth', 2);
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    title('x', 'FontSize', 10);
    xlabel('Frequency (Hz)', 'FontSize', 10);
    ylabel('Power (a.u.)', 'FontSize', 10);
end
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% For each channel...
for i = 1:SIGMA.nbchan+1
    if i <= SIGMA.nbchan
        meas = MEAS(idx_freqfit, i);
        chanlbl = lower(SIGMA.chanlocs(i).labels);
    else
        meas = mean(MEAS(idx_freqfit, :), 2);
        chanlbl = 'avg';
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Fit gaussian curves, or use 'smoothdata'
    switch fitMethod
        case 'gaussianfit'
            freq = ascolumn(FREQ(idx_freqfit));
            fitpow = ascolumn(meas);
            try
            [fitresult, offset] = gaussianfit(freq, fitpow, ngaussians);
            catch ME
                keyboard
            end
            fitpow = fitresult(freq) + offset;
            span = nan;
            coeffvals = coeffvalues(fitresult);
            mu = coeffvals(2);
            [~, idx_mu] = min(abs(freq - mu));
            peak = fitpow(idx_mu);
            bandwidth = coeffvals(3);
            lolim = mu - bandwidth/2;
            hilim = mu + bandwidth/2;
        case 'smoothdata'
            span = 1;
            freq = ascolumn(FREQ(idx_freqfit));
            fitpow = smoothdata(meas, 'gaussian', span);
            pks = findpeaks(fitpow);
            while length(pks) > 1 && span < length(meas)
                span = span + 1;
                fitpow = smoothdata(meas, 'gaussian', span);
                pks = findpeaks(fitpow(fitpow > mean(fitpow)));
            end
            fitpow(freq < 0.0075) = 0;
            [peak, idx_max] = max(fitpow);
            mu = freq(idx_max);
            idx_edge = idx_max;
            while fitpow(idx_edge) > fitpow(idx_edge+1)
                idx_edge = idx_edge+1;
                if idx_edge == length(fitpow)
                    break
                end
            end
            lolim = mu - (freq(idx_edge) - mu);
            hilim = freq(idx_edge);
            bandwidth = hilim - lolim;
            coeffvals = NaN;
            offset = 0;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Store spectral data and infraslow modulation peak
    FIT(i).freq = FREQ(idx_freqfit);
    FIT(i).pow = meas;
    FIT(i).fitpow = fitpow;
    FIT(i).mu = mu;
    FIT(i).peak = peak;
    FIT(i).lower = lolim;
    FIT(i).upper = hilim;
    FIT(i).bandwidth = bandwidth;
    FIT(i).coeffvals = coeffvals;
    % ---------------------------------------------------------------------
    % Create figure
    if doPlot
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        set(H.data, 'XData', FIT(i).freq, 'YData', FIT(i).pow)
        switch fitMethod
            case 'gaussianfit'
                set(H.span, 'XData', [NaN, NaN])
                for j = 1:length(H.gaus)
                    H.gaus(j).XData = NaN;
                    H.gaus(j).YData = NaN;
                end
                for j = 1:3:length(coeffvals)
                    a = coeffvals(j);
                    b = coeffvals(j+1);
                    c = coeffvals(j+2);
                    H.gaus((2+j)/3).XData = freq;
                    H.gaus((2+j)/3).YData = a.*exp(-1 .* ((freq - b).^2)/ (2.*c.^2));
                end
            case 'smoothdata'
                set(H.span, 'XData', [FIT(i).freq(end-span), FIT(i).freq(end)])
        end
        set(H.fit, 'XData', FIT(i).freq, 'YData', FIT(i).fitpow)
        set(H.thres, 'XData', [0.0075, 0.0075], 'YData', [0, max(meas(:))*1.2])
        set(H.lolim, 'XData', [FIT(i).lower, FIT(i).lower], 'YData', [0, max(meas(:))*1.2])
        set(H.hilim, 'XData', [FIT(i).upper, FIT(i).upper], 'YData', [0, max(meas(:))*1.2])
        set(H.peak, 'XData', FIT(i).mu, 'YData', FIT(i).peak)
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        Ax.Title.String = sprintf('Channel %s', chanlbl);
        Ax.XLim = [0, 0.1];
        Ax.YLim = [0, FIT(i).peak*1.2];
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        kv = filename2struct(SIGMA.setname);
        if contains(kv.desc, 'spindle')
            kv.desc = sprintf('desc-infraslowspindle%s', chanlbl);
            outfilepath = [SIGMA.filepath, '/images/infraslowspindle'];
        else
            kv.desc = sprintf('desc-infraslowsigma%s', chanlbl);
            outfilepath = [SIGMA.filepath, '/images/infraslow'];
        end
        outfilename = [struct2filename(kv), '.png'];
        if exist(outfilepath, 'dir') == 0
            mkdir(outfilepath);
        end
        exportgraphics(Fig, fullfile(outfilepath, outfilename));
    end
end
% -------------------------------------------------------------------------
% Create topoplots
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Colormap
if doPlot
    CMap = load('colormap_roma.mat');
    for fld = {'peak', 'mu'}
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        % YData
        YData = [FIT.(fld{:})];
        YData(end) = []; % Remove average
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        % Create figure
        Fig = figure('Color', 'w', 'Position', [560, 640, 290 200]);
        Ax = axes();
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        % Draw topoplot
        switch fld{:}
            case 'peak'
                contourvals = (YData > 1) + (YData > 2);
                CLim = [0, max(YData(:))];
            case 'mu'
                contourvals = (YData > 0.01) + (YData > 0.02) + (YData > 0.03);
                CLim = [0.0075, 0.04];
        end
        topoplot(YData, SIGMA.chanlocs, ...
            'headrad', 0.575, ...
            'whitebk', 'on', ...
            'conv', 'on', ...
            'numcontour', length(unique(contourvals))-1, ...
            'contourvals', contourvals);
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        % Reset figure and axis properties
        Fig.Color = 'w';
        Fig.Colormap = CMap.roma;
        Ax.CLim = CLim;
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        % Colobar
        CBar = colorbar();
        switch fld{:}
            case 'peak'
                CBar.Ticks = sort([0, 1, 2, max(YData(:))]);
                CBar.Label.String = 'peak amplitude (norm.)';
            case 'mu'
                CBar.Ticks = sort([0, 0.01, 0.02, 0.03, max(YData(:))]);
                CBar.Label.String = 'mean frequency (Hz)';
        end
        CBar.FontSize = 8;
        CBar.Label.FontSize = 10;
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        kv = filename2struct(SIGMA.setname);
        if contains(kv.desc, 'spindle')
            kv.desc = sprintf('desc-topospindle%s', fld{:});
            outfilepath = [SIGMA.filepath, '/images/infraslowspindle'];
        else
            kv.desc = sprintf('desc-topowsigma%s', fld{:});
            outfilepath = [SIGMA.filepath, '/images/infraslow'];
        end
        outfilename = [struct2filename(kv), '.png'];
        if exist(outfilepath, 'dir') == 0
            mkdir(outfilepath);
        end
        exportgraphics(Fig, fullfile(outfilepath, outfilename));
    end
end
close all

