function css_crosscorr(EEG, cfg)

% Load the visual inspection file
excl = readtable('inspect/ecg/hr_excl_nrembout.xlsx');
maxlag = 120;
% Find accompanying HR dataset
kv = filename2struct(EEG.setname);
hrfilename = dir(sprintf('derivatives/EEG-segmented/sub-%s/ses-%s/sub-%s_ses-%s_task-psg_desc-hrnrembout_hr.set', kv.sub, kv.ses, kv.sub, kv.ses));
ecgfilename = dir(sprintf('derivatives/EEG-segmented/sub-%s/ses-%s/sub-%s_ses-%s_task-psg_desc-ecgnrembout_ecg.set', kv.sub, kv.ses, kv.sub, kv.ses));
if isempty(hrfilename)
    return
end
HR = LoadDataset(fullfile(hrfilename(1).folder, hrfilename(1).name), 'all');
HR = pop_resample(HR, EEG.srate);
HR.times = linspace(HR.xmin, HR.xmax, HR.pnts);

ECG = LoadDataset(fullfile(ecgfilename(1).folder, ecgfilename(1).name), 'all');
ECG.times = linspace(ECG.xmin, ECG.xmax, ECG.pnts);

% Init output structure
X = struct();
X.data = [];
% Get the sample indices of the NREM bouts
bouts = find(strcmpi({EEG.event.type}, 'boundary'));
bouts = [1, ceil([EEG.event(bouts).latency]), EEG.pnts];

% Run across all channels
tr = now();
for c = 1:EEG.nbchan
    % Initialize empty R matrix
    R = [];
    close all
    figure
    hold on
    for b = 1:length(bouts)-1
        % Check if we need to skip this bout
        if any(strcmpi(kv.sub, excl.sub) & strcmpi(kv.ses, excl.ses) & (excl.excl_bout == b | excl.excl_bout == 999))
            continue
        end
        % Set the sample indices for this bout
        this_bout = [bouts(b), bouts(b+1)];
        % Generate windows of 'maxlag' seconds
        wins = ascolumn(this_bout(1):6*HR.srate:this_bout(2)-maxlag*HR.srate);
        wins = [wins, wins+maxlag*HR.srate-1]; %#ok<AGROW> 
        for w = 1:size(wins, 1)
            % Get the HR and Sigma power timeseries data of this bout
            x = nanzscore(HR.data(1, wins(w, 1):wins(w, 2)));
            y = nanzscore(EEG.data(c, wins(w, 1):wins(w, 2)));
            % Remove NaN's
            idx_rm = isnan(x) | isnan(y);
            x(idx_rm) = mean(x, 'omitnan');
            y(idx_rm) = mean(y, 'omitnan');
            if isempty(x)
                keyboard
                error('No HR data')
            end
            % Calculate cross corr and store in matrix
            [r, lags] = xcorr(detrend(x), detrend(y), maxlag*HR.srate/2, 'Normalized');
            R = [R; r]; %#ok<AGROW> 
        end
        tr = remainingTime(tr, EEG.nbchan*length(bouts)-1);
    end
    % if there is no data in 'R' then there were no good quality bouts and
    % we cannot get the cross correlations, stop here
    if isempty(R)
        return
    end
    X.data = [X.data; mean(R, 'omitnan')];
end

% Set times vector
X.times = lags./EEG.srate;
X.xmin = X.times(1);
X.xmax = X.times(end);
X.pnts = length(X.times);
X.nbchan = size(X.data, 1);
% Extract featres
Features = struct();
% Get max within +/- 50 seconds
idx_lag = X.times > -50 & X.times < 50;
[xcorr_peak_amp, idx_max] = max(X.data(:, idx_lag)'); %#ok<UDIM>
idx_max = idx_max+find(X.times > -50, 1, 'first')-1;
% Lag
Features(1).label = 'xcorr_lag';
Features(1).type = 's';
Features(1).data = X.times(idx_max);
% Peak amplitude
Features(2).label = 'xcorr_peak';
Features(2).type = 'r';
Features(2).data = xcorr_peak_amp;
% FWHM
fwhm = nan(size(X.data, 1), 1);
for c = 1:size(X.data, 1)
    % Find the leading edge
    lft = fliplr(X.data(c, 1:idx_max(c)));
    lft = find(lft < xcorr_peak_amp(c)/2, 1, 'first');
    % Find the falling edge
    rgt = X.data(c, idx_max(c):end);
    rgt = find(rgt < xcorr_peak_amp(c)/2, 1, 'first');
    fwhm(c) = lft+rgt-1;
end
Features(3).label = 'xcorr_fwhm';
Features(3).type = 's';
Features(3).data = fwhm./HR.srate;

% Save output
css_createfstlvloutput(cfg.outfilepath, Features, X);

% Plotting (if requested)
if cfg.doPlot
    close all
    Fig = figure('Color', 'w', 'Position', [1630, 775, 290 200]);
    Ax = axes('NextPlot', 'add', ...
        'FontSize', 10, ...
        'XTick', [X.times(1), 0, X.times(end)], ...
        'YLim', [-1.*max(abs(X.data(:))), 1.*max(abs(X.data(:)))].*1.1, ...
        'YTick', -1:0.1:1, ...
        'XGrid', 'on', ...
        'YGrid', 'on', ...
        'TickLength', [0 0], ...
        'Box', 'on', ...
        'Color', [0.94 0.95 0.97]);
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    H.data = plot(Ax, NaN, NaN, '-', 'Color', [0.6 0.6 0.6]);
    H.peak = plot(Ax, NaN, NaN, '.r', 'MarkerSize', 8);
    H.fwhm = plot(Ax, NaN, NaN, '-r');
    for c = 1:size(X.data, 1)
        H.data.XData = X.times;
        H.data.YData = X.data(c, :);
        H.peak.XData = Features(1).data(c);
        H.peak.YData = Features(2).data(c);
        H.fwhm.XData = [Features(1).data(c) - Features(3).data(c)./2, Features(1).data(c) + Features(3).data(c)./2];
        H.fwhm.YData = [Features(2).data(c), Features(2).data(c)]./2;
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        kv = filename2struct(reverse_fileparts(cfg.outfilepath));
        kv.desc = [kv.desc, 'e', num2str(c)];
        kv.filetype = 'xcorrinspect.png';
        figoutfilepath = fullfile(fileparts(cfg.outfilepath), 'images', struct2filename(kv));
        if exist([fileparts(cfg.outfilepath), '/images'], 'dir') == 0
            mkdir([fileparts(cfg.outfilepath), '/images']);
        end
        exportgraphics(Fig, figoutfilepath);
    end
end
end