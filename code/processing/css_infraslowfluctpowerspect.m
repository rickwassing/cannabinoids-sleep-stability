function [FIT] = css_infraslowfluctpowerspect(EEG, cfg)
% -------------------------------------------------------------------------
% Constants
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Upper limit of infraslow modulation frequency:
freq_upper_limit = 0.5;
% What method to use to extract ISF peaks?
fitMethod = 'gaussianfit'; %  Options 'gaussianfit' or 'smoothdata'
% -------------------------------------------------------------------------
% Extract bouts
idx_bouts = round([1, EEG.event(strcmpi({EEG.event.type}, 'boundary')).latency, EEG.pnts+1]);
bouts = {};
for i = 1:length(idx_bouts)-1
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Extract the bout and demean the data
    this_idx = idx_bouts(i):(idx_bouts(i+1)-1);
    x = EEG.data(:, this_idx)';
    % Remove nan's
    x(sum(isnan(x), 2) > 0, :) = [];
    % High-pass
    if cfg.doFilter
        forder = pop_firwsord('hamming', EEG.srate, 0.01);
        if length(x) <= forder*3
            x = detrend(x);
        else
            b = fir1(forder, 0.005/(EEG.srate/2), 'high');
            x = filtfilt(b, 1, double(x));
        end
    else
        x = detrend(x, 0); % demean
    end
    bouts{i} = x; %#ok<AGROW>
end
% -------------------------------------------------------------------------
% Init output
nfft = 300.*EEG.srate; % 2*N-1
window = hanning(nfft);
overlap = round(nfft/2);
FREQS = linspace(0, EEG.srate/2, nfft/2+1);
idx_f = find(FREQS <= freq_upper_limit);
% -------------------------------------------------------------------------
% First and last bouts may have NaN values because of the edge artefacts of
% the wavelet analysis (half the wavelet length). Remove bouts that are too
% short.
idx_rm = cellfun(@(b) length(b), bouts) < nfft;
bouts(idx_rm) = [];
ABSPXX = zeros(length(bouts), length(idx_f), EEG.nbchan);
% -------------------------------------------------------------------------
% For each bout in the EEG data...
for i = 1:length(bouts)
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Apply PWELCH to the freq-band power timeseries
    [pwel_infraslow, pwel_freqs] = pwelch(bouts{i}, window, overlap, nfft, EEG.srate, 'power'); %#ok<ASGLU>
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Store in matrix
    ABSPXX(i, :, :) = pwel_infraslow(idx_f, :, :);
end
% -------------------------------------------------------------------------
% Crop frequency range
FREQS = FREQS(idx_f);
% -------------------------------------------------------------------------
% Calculate relative power (normalize to total power below 0.1 Hz)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
ABSPXX = squeeze(mean(ABSPXX, 1)); % Average across bouts
grand_total_power = mean(mean(abs(ABSPXX(FREQS > 0.0075 & FREQS <= 0.1, :)), 1));
RELPXX = ABSPXX ./ grand_total_power;
% -------------------------------------------------------------------------
% We're fitting curves to frequencies under 0.1 Hz
idx_freqfit = FREQS <= 0.1;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Init the output structure
FIT = struct();
% Init figure
if cfg.doPlot
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    close all
    Fig = figure('Color', 'w', 'Position', [560, 640, 290 200]);
    Ax = axes('NextPlot', 'add', ...
        'FontSize', 10, ...
        'YTick', 1:10, ...
        'XGrid', 'on', ...
        'TickLength', [0 0], ...
        'Box', 'on', ...
        'Color', [0.94 0.95 0.97]);
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
for i = 1:EEG.nbchan
    if i <= EEG.nbchan
        abspxx = ABSPXX(idx_freqfit, i);
        relpxx = RELPXX(idx_freqfit, i);
        chanlbl = lower(EEG.chanlocs(i).labels);
    else
        abspxx = mean(ABSPXX(idx_freqfit, :), 2);
        relpxx = mean(RELPXX(idx_freqfit, :), 2);
        chanlbl = 'avg';
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Fit gaussian curves, or use 'smoothdata'
    if cfg.normalize
        this_fit = fitisfspect(fitMethod, ascolumn(FREQS(idx_freqfit)), ascolumn(relpxx));
    else
        this_fit = fitisfspect(fitMethod, ascolumn(FREQS(idx_freqfit)), ascolumn(abspxx));
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Store spectral data and infraslow modulation peak
    FIT(i).freq = this_fit.freq;
    FIT(i).pow = this_fit.pow;
    FIT(i).fitpow = this_fit.fitpow;
    FIT(i).mu = this_fit.mu;
    FIT(i).peak = this_fit.peak;
    FIT(i).lower = this_fit.lower;
    FIT(i).upper = this_fit.upper;
    FIT(i).bandwidth = this_fit.bandwidth;
    FIT(i).coeffvals = this_fit.coeffvals;
    FIT(i).span = this_fit.span;
    FIT(i).offset = this_fit.offset;
    % ---------------------------------------------------------------------
    % Create figure
    if cfg.doPlot
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        set(H.data, 'XData', FIT(i).freq, 'YData', FIT(i).pow)
        switch fitMethod
            case 'gaussianfit'
                set(H.span, 'XData', [NaN, NaN])
                for j = 1:length(H.gaus)
                    H.gaus(j).XData = NaN;
                    H.gaus(j).YData = NaN;
                end
                for j = 1:3:length(FIT(i).coeffvals)
                    a = FIT(i).coeffvals(j);
                    b = FIT(i).coeffvals(j+1);
                    c = FIT(i).coeffvals(j+2);
                    H.gaus((2+j)/3).XData = FIT(i).freq;
                    H.gaus((2+j)/3).YData = a.*exp(-1 .* ((FIT(i).freq - b).^2)/ (2.*c.^2));
                end
            case 'smoothdata'
                set(H.span, 'XData', [FIT(i).freq(end - FIT(i).span), FIT(i).freq(end)])
        end
        set(H.fit, 'XData', FIT(i).freq, 'YData', FIT(i).fitpow)
        set(H.lolim, 'XData', [FIT(i).lower, FIT(i).lower], 'YData', [0, FIT(i).peak])
        set(H.hilim, 'XData', [FIT(i).upper, FIT(i).upper], 'YData', [0, FIT(i).peak])
        set(H.peak, 'XData', FIT(i).mu, 'YData', FIT(i).peak)
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        Ax.Title.String = sprintf('Channel %s', chanlbl);
        Ax.XLim = [FIT(i).freq(2), max(FIT(i).freq)];
        Ax.YLim = [0, FIT(i).peak*1.2];
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        kv = filename2struct(reverse_fileparts(cfg.outfilepath));
        kv.desc = [kv.desc, chanlbl];
        kv.filetype = 'isfspect.png';
        figoutfilepath = fullfile(fileparts(cfg.outfilepath), 'images', struct2filename(kv));
        if exist([fileparts(cfg.outfilepath), '/images'], 'dir') == 0
            mkdir([fileparts(cfg.outfilepath), '/images']);
        end
        exportgraphics(Fig, figoutfilepath);
    end
end
% -------------------------------------------------------------------------
% Create topoplots
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Colormap
if cfg.doPlot
    CMap = load('colormap_roma.mat');
    for fld = {'peak', 'mu', 'bandwidth'}
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        % YData
        YData = [FIT.(fld{:})];
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
            case 'bandwidth'
                contourvals = (YData > 0.005) + (YData > 0.010) + (YData > 0.015) + (YData > 0.020) + (YData > 0.025);
                CLim = [0, max(YData(:))];
        end
        topoplot(YData, EEG.chanlocs, ...
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
                CBar.Ticks = unique([0, 1, 2, max(YData(:))]);
                CBar.Label.String = 'peak amplitude (norm.)';
            case 'mu'
                CBar.Ticks = unique([0, 0.01, 0.02, 0.03, max(YData(:))]);
                CBar.Label.String = 'mean frequency (Hz)';
            case 'bandwidth'
                CBar.Ticks = unique([0, 0.005, 0.010, 0.015, 0.020, 0.025, max(YData(:))]);
                CBar.Label.String = 'bandwidth (Hz)';
        end
        CBar.FontSize = 8;
        CBar.Label.FontSize = 10;
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        kv = filename2struct(reverse_fileparts(cfg.outfilepath));
        kv.desc = [kv.desc, fld{:}];
        kv.filetype = 'isftopo.png';
        figoutfilepath = fullfile(fileparts(cfg.outfilepath), 'images', struct2filename(kv));
        if exist([fileparts(cfg.outfilepath), '/images'], 'dir') == 0
            mkdir([fileparts(cfg.outfilepath), '/images']);
        end
        exportgraphics(Fig, figoutfilepath);
    end
end
close all
% -------------------------------------------------------------------------
% Initialize first level output
ISF = struct();
ISF.data = cat(2, FIT.pow)';
ISF.freqs = FIT(1).freq;
ISF.freqstep = mean(diff(FIT(1).freq));
% Extract the features used in group level analysis
ISF.features = css_extractfeatures(FIT, ISF);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
css_createfstlvloutput(strrep(cfg.outfilepath, '.set', '.mat'), ISF.features, ISF);

end