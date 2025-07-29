% -------------------------------------------------------------------------
% Load all REM onsets and plot the ~100 seconds of sigma ISF
SigmaFiles = dir('derivatives/EEG-segmented/sub-*/ses-placebo/sub-*-sigmaprerembout_pow.set');
% -------------------------------------------------------------------------
% Settings
chanlocs = template_to_chanlocs(which('GSN-HydroCel-257.sfp'));
[~, system] = ishdeeg({chanlocs.labels});
incl = hdeeg_scalpchannels(system);
chanlocs = chanlocs(ismember({chanlocs.labels}, incl));
chanlocs = channel_clusters(chanlocs, 'mff');
idx_chan = [77, 78, 79, 85, 86, 87, 92, 93, 94, 99, 104, 110, 111, 112, 113, 120, 121, 122];
chanlabels = {chanlocs(idx_chan).labels};
append_type = 'zeros';
ford_mult = 2;
% -------------------------------------------------------------------------
% Set filter config
SIGMA = LoadDataset(fullfile(SigmaFiles(1).folder, SigmaFiles(1).name), 'all');
filtcfg.cutoff = [0.0138, 0.0259];
filtcfg.wintype = 'kaiser';
filtcfg.transbw = 1/50;
filtcfg.rippledev = 0.05;
filtcfg.warg = 2;
filtcfg.adj = [0 0];
filtcfg.order = pop_firwsord(filtcfg.wintype, SIGMA(1).srate, filtcfg.transbw, filtcfg.rippledev); % transition bandwidth of 0.01 Hz
filtcfg.order
% -------------------------------------------------------------------------
% Settings for detecting troughs and peaks in the sigma ISF
evcfg = struct();
evcfg.srate = SIGMA(1).srate;
evcfg.pkdist = 30;
evcfg.pkheight = 0.75.*pi;
evcfg.mindelay = 20; % seconds. Minimum time prior to REM onset to identify the ISF trough
% -------------------------------------------------------------------------
% Init
H = struct();
cnt = 0;
FName = [];
XData = [];
YData = [];
% -------------------------------------------------------------------------
% Load and process files
for i = 2:length(SigmaFiles)
    % ---------------------------------------------------------------------
    % Load file and select parietal channels only
    SIGMA = LoadDataset(fullfile(SigmaFiles(i).folder, SigmaFiles(i).name), 'all');
    SIGMA = pop_select(SIGMA, 'channel', idx_chan);
    % ---------------------------------------------------------------------
    % Extract bout indexes
    bouts = ([SIGMA.event(strcmpi({SIGMA.event.type}, 'boundary')).latency]);
    bouts = [0.5, bouts, SIGMA.pnts+0.5];
    % For each bout
    for b = 1:length(bouts)-1
        % -----------------------------------------------------------------
        % Cut and append the data
        sigma = pop_select(SIGMA, 'point', [bouts(b), bouts(b+1)]);
        sigma.data = zscoreacrosschannels(sigma.data);
        sigma.data = detrend(sigma.data', 0, 'omitnan')'; % demean
        sigma = signalappend(sigma, append_type, ford_mult, filtcfg);
        % -----------------------------------------------------------------
        % Prepend and append data to cover for filter edge artefact
        sigma = executeappending(sigma);
        % -----------------------------------------------------------------
        % Filter the data
        sigma = pop_firws(sigma, ...
            'fcutoff', filtcfg.cutoff+filtcfg.adj, ...
            'ftype', 'bandpass', ...
            'wtype', filtcfg.wintype, ...
            'warg', filtcfg.warg, ...
            'forder', filtcfg.order, ...
            'plotfresp', false, ...
            'minphase', 0);
        % -----------------------------------------------------------------
        % Remove appending
        sigma = executeappending(sigma, 'remove');
        % -----------------------------------------------------------------
        % Adjust time vector where zero indicates REM onset
        tzero = sigma.event(find(strcmpi({sigma.event.type}, 'remeps'), 1, 'first')).latency;
        remepdur = sigma.event(find(strcmpi({sigma.event.type}, 'remeps'), 1, 'first')).duration;
        sigma.times = ((0:sigma.pnts-1)-tzero)./sigma.srate;
        % -----------------------------------------------------------------
        % Apply Hilbert
        cnt = cnt+1;
        H(cnt).sigma = sigma.data';
        H(cnt).sigma = detrend(H(cnt).sigma, 0);
        H(cnt).times = sigma.times;
        H(cnt).cmx = hilbert(H(cnt).sigma);
        H(cnt).ang = angle(H(cnt).cmx);
        H(cnt).amp = abs(H(cnt).cmx);
        H(cnt).remepdur = remepdur;
        % -----------------------------------------------------------------
        % Find peaks and troughs at the REM transition, and the ones before that
        H(cnt).event = prerempeakstroughs(H(cnt), evcfg);
    end
end
disp('done')
%%

close all;
Fig = figure('Position', [680, 87, 777, 890]);
Ax = axes(Fig, 'NextPlot', 'add');

a = 10;
b = 10;
m = 20;

add = 30;
n = 15;
chan = 2;


for i = (1:n)+add
    % Plot timeseries and phase angle
    idx_tr = H(i).event.idx_tr{H(i).event.channel == chan & strcmpi(H(i).event.stage, 'rem')};
    plot(H(i).times(idx_tr), H(i).sigma(idx_tr, chan)*m+i*b*pi+a, 'o', 'Color', [0.6 0.4 0.9])
    plot(H(i).times, H(i).sigma(:, chan)*m+i*b*pi+a, '-', 'Color', [0.6 0.4 0.9])
    plot(H(i).times, H(i).ang(:, chan)+i*b*pi, '-k')
    % Plot markers for the ISF prior to REM transitions
    idx_ev = H(i).event.channel == chan & strcmpi(H(i).event.stage, 'nrem');
    idx_tr = H(i).event.latency(idx_ev);
    idx_pk = H(i).event.peak_latency(idx_ev);
    plot(H(i).times(idx_tr), H(i).sigma(idx_tr, chan)*m+i*b*pi+a, '.b', 'MarkerSize', 6)
    plot(H(i).times(idx_pk), H(i).sigma(idx_pk, chan)*m+i*b*pi+a, '.b', 'MarkerSize', 6)
    for j = 1:length(idx_tr)
        plot([H(i).times(idx_tr(j)), H(i).times(idx_tr(j)), H(i).times(idx_pk(j))], [H(i).sigma(idx_tr(j), chan)*m+i*b*pi+a, H(i).sigma(idx_pk(j), chan)*m+i*b*pi+a, H(i).sigma(idx_pk(j), chan)*m+i*b*pi+a], '-b')
    end
    % Plot markers for the ISF at REM transitions
    idx_ev = H(i).event.channel == chan & strcmpi(H(i).event.stage, 'rem');
    idx_tr = H(i).event.latency(idx_ev);
    idx_pk = H(i).event.peak_latency(idx_ev);
    plot(H(i).times(idx_tr), H(i).sigma(idx_tr, chan)*m+i*b*pi+a, '.r', 'MarkerSize', 6)
    plot(H(i).times(idx_pk), H(i).sigma(idx_pk, chan)*m+i*b*pi+a, '.r', 'MarkerSize', 6)
    plot([H(i).times(idx_tr), H(i).times(idx_tr), H(i).times(idx_pk)], [H(i).sigma(idx_tr, chan)*m+i*b*pi+a, H(i).sigma(idx_pk, chan)*m+i*b*pi+a, H(i).sigma(idx_pk, chan)*m+i*b*pi+a], '-r')

    if i == length(H)
        break
    end
end

Ax.TickLength = [0 0];
Ax.YTick = ((1:n)+add).*b.*pi;
Ax.YTickLabel = (1:n)+add;
Ax.YLim = [Ax.YTick(1), Ax.YTick(end)]+[-2.*pi 2.*pi+a];
Ax.XTick = -300:50:100;

plot([-evcfg.mindelay, -evcfg.mindelay], [Ax.YTick(1), Ax.YTick(end)]+[-2.*pi 2.*pi+a], ':r')
plot([0, 0], [Ax.YTick(1), Ax.YTick(end)]+[-2.*pi 2.*pi+a], ':k')
plot([30, 30], [Ax.YTick(1), Ax.YTick(end)]+[-2.*pi 2.*pi+a], ':k')

%%
close all

XData_nrem = arrayfun(@(h) h.event.latency(strcmpi(h.event.stage, 'nrem')), H, 'UniformOutput', false);
XData_rem = arrayfun(@(h) h.event.latency(strcmpi(h.event.stage, 'rem')), H, 'UniformOutput', false);
Dur_nrem = arrayfun(@(h) h.event.duration(strcmpi(h.event.stage, 'nrem')), H, 'UniformOutput', false);
Dur_rem = arrayfun(@(h) h.event.duration(strcmpi(h.event.stage, 'rem')), H, 'UniformOutput', false);
Amp_nrem = arrayfun(@(h) h.event.amplitude(strcmpi(h.event.stage, 'nrem')), H, 'UniformOutput', false);
Amp_rem = arrayfun(@(h) h.event.amplitude(strcmpi(h.event.stage, 'rem')), H, 'UniformOutput', false);

Av_XData_nrem = arrayfun(@(h) mean(h.event.latency(strcmpi(h.event.stage, 'nrem'))), H, 'UniformOutput', true);
Av_XData_rem = arrayfun(@(h) mean(h.event.latency(strcmpi(h.event.stage, 'rem'))), H, 'UniformOutput', true);
Av_Dur_nrem = arrayfun(@(h) mean(h.event.duration(strcmpi(h.event.stage, 'nrem'))), H, 'UniformOutput', true);
Av_Dur_rem = arrayfun(@(h) mean(h.event.duration(strcmpi(h.event.stage, 'rem'))), H, 'UniformOutput', true);
Av_Amp_nrem = arrayfun(@(h) mean(h.event.amplitude(strcmpi(h.event.stage, 'nrem'))), H, 'UniformOutput', true);
Av_Amp_rem = arrayfun(@(h) mean(h.event.amplitude(strcmpi(h.event.stage, 'rem'))), H, 'UniformOutput', true);

Dur_remep = [H.remepdur];

Fig = figure();
Ax = axes(Fig, 'NextPlot','add');
plot(cat(1, XData_nrem{:}), cat(1, Dur_nrem{:}), '.b')
plot(cat(1, XData_rem{:}), cat(1, Dur_rem{:}), '.r')
Ax.Title.String = 'duration';

Fig = figure();
Ax = axes(Fig, 'NextPlot','add');
plot(cat(1, XData_nrem{:}), cat(1, Amp_nrem{:}), '.b')
plot(cat(1, XData_rem{:}), cat(1, Amp_rem{:}), '.r')
Ax.Title.String = 'amplitude';

Fig = figure();
Ax = axes(Fig, 'NextPlot','add');
histogram(cat(1, Dur_nrem{:}), 'Normalization', 'probability', 'BinEdges', 0:20:1000)
histogram(cat(1, Dur_rem{:}), 'Normalization', 'probability', 'BinEdges', 0:20:1000)
Ax.Title.String = 'duration';

Fig = figure();
Ax = axes(Fig, 'NextPlot','add');
histogram(cat(1, Amp_nrem{:}), 50, 'Normalization', 'probability', 'BinEdges', 0:0.05:2)
histogram(cat(1, Amp_rem{:}), 50, 'Normalization', 'probability', 'BinEdges', 0:0.05:2)
Ax.Title.String = 'amplitude';


Fig = figure();
Ax = axes(Fig, 'NextPlot','add');
plot(Av_XData_nrem, Av_Dur_nrem, '.b')
plot(Av_XData_rem, Av_Dur_rem, '.r')
Ax.Title.String = 'average duration';

Fig = figure();
Ax = axes(Fig, 'NextPlot','add');
plot(Av_XData_nrem, Av_Amp_nrem, '.b')
plot(Av_XData_rem, Av_Amp_rem, '.r')
Ax.Title.String = 'average amplitude';

Fig = figure();
Ax = axes(Fig, 'NextPlot','add');
histogram(Av_Dur_nrem, 'Normalization', 'probability', 'BinEdges', 0:20:1000)
histogram(Av_Dur_rem, 'Normalization', 'probability', 'BinEdges', 0:20:1000)
[~, pval, ~, stats] = ttest(Av_Dur_nrem, Av_Dur_rem);
Ax.Title.String = sprintf('average duration t = %.2f, p = %.2f', stats.tstat, pval);

Fig = figure();
Ax = axes(Fig, 'NextPlot','add');
histogram(Av_Amp_nrem, 'Normalization', 'probability', 'BinEdges', 0:0.05:2)
histogram(Av_Amp_rem, 'Normalization', 'probability', 'BinEdges', 0:0.05:2)
[~, pval, ~, stats] = ttest(Av_Amp_nrem, Av_Amp_rem);
Ax.Title.String = sprintf('average amplitude t = %.2f, p = %.2f', stats.tstat, pval);

Fig = figure();
Ax = axes(Fig, 'NextPlot','add');
plot(Av_Dur_rem, Dur_remep, '.r')
lsline
Ax.XLabel.String = 'ISF Duration';
Ax.YLabel.String = 'REM Duration';
