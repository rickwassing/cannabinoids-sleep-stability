function css_plot_1c(type, varargin)

% Load group-level results, colormap, and chanlocs
load(sprintf('group-level/a1c_pairttest_cmass/%ssigma/glm.mat', type))
grpx = load('group-level/a1c_pairttest_cmass/xcorr/glm.mat');
GLMXC = grpx.GLM;
load('colormap_roma.mat')
load('colormap_batlow.mat')
BOUTS.a = readtable('group-level/nrembout_number.csv');
BOUTS.b = readtable('group-level/nrembout_duration.csv');
chanlocs = template_to_chanlocs(which('GSN-HydroCel-257.sfp'));
[~, system] = ishdeeg({chanlocs.labels});
incl = hdeeg_scalpchannels(system);
chanlocs = chanlocs(ismember({chanlocs.labels}, incl));
chanlocs = channel_clusters(chanlocs, 'mff');
% Set limits
switch type
    case 'abs'
        AmpYLim = [0, 2];
    case 'norm'
        AmpYLim = [0, 4];
end
if nargin > 1
    idx_chan = find(ismember({chanlocs.labels}, varargin{1}));
else
    idx_chan = 1:length(chanlocs);
end
% Load example dataset
SIGMA = LoadDataset('derivatives/EEG-segmented/sub-r005/ses-etc120/sub-r005_ses-etc120_task-psg_desc-sigmanrembout_pow.set', 'all');
EEG = LoadDataset('derivatives/EEG-processed/sub-r005/ses-etc120/sub-r005_ses-etc120_task-psg_desc-sigma_pow.set', 'all');
HR = LoadDataset('derivatives/EEG-preproc/sub-r005/ses-etc120/sub-r005_ses-etc120_task-psg_desc-preprochr_hr.set', 'all');
bouts = getnrembouts(css_eeglab2hypnogram(EEG), EEG.srate, 300);
EEG = pop_select(EEG, 'time', [0, bouts(end, 2)+300]);
EEG.times = linspace(EEG.xmin, EEG.xmax, EEG.pnts);
HR = pop_select(HR, 'time', [0, bouts(end, 2)+300]);
HR.times = linspace(HR.xmin, HR.xmax, HR.pnts);
% Load data from both groups
Files = dir(sprintf('derivatives/EEG-output-fstlvl/sub-*/ses-*/sub-*%ssigma*interp_fstlvl.mat', type));
ISF = [];
for i = 1:length(Files)
    if i == 1
        ISF = LoadDataset(fullfile(Files(i).folder, Files(i).name), 'matrix');
    else
        ISF(i) = LoadDataset(fullfile(Files(i).folder, Files(i).name), 'matrix');
    end
end
% Cross correlation between Sigma power and HR timeseries
Files = dir('derivatives/EEG-output-fstlvl/sub-*/ses-*/sub-*_desc-nremboutxcorr120s_fstlvl.mat');

XC = [];
for i = 1:length(Files)
    tmp = LoadDataset(fullfile(Files(i).folder, Files(i).name), 'matrix');
    if isempty(XC)
        XC = tmp;
    else
        XC(i) = tmp;
    end
end
% -------------------------------------------------------------------------
% Create new figure
close all
clear Ax
ai = 0;
Fig = figure('Color', 'w');
Fig.Units = 'centimeters';
Fig.Position = [5, 15, 18, 8];
% -------------------------------------------------------------------------
% Example hypnogram
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Create axes
ai = ai+1;
Ax(ai) = axes(...
    'NextPlot', 'add', ...
    'FontSize', 10, ...
    'Box', 'on', ...
    'OuterPosition', [0, 1-3/12, 8/12, 3/12] + [-0.03 0.035 0.03 -0.035]);
plotHypnogram(Ax(ai), EEG, 'DoPlotCycles', false);
selbout = 4;
for i = 1:size(bouts)
    YData = [-3.75; -3.75; 1.75; 1.75; -3.75];
    XData = [bouts(i, 1); bouts(i, 2); bouts(i, 2); bouts(i, 1); bouts(i, 1)] ./ (24*60*60);
    Vertices = [XData, YData];
    Faces = 1:4;
    patch(Ax, 'Faces', Faces, 'Vertices', Vertices, ...
        'EdgeColor', [0, 0, 0], ...
        'FaceColor', standard_colors('blue'), ...
        'FaceAlpha', 0.15, ...
        'EdgeAlpha', 0.30, ...
        'LineStyle', ':');
end
% Zoom lines
plot(Ax(ai), [bouts(selbout, 1), 0]./(24*60*60), [-3.9 -6], ':', 'Color', [0.8, 0.8, 0.8], 'LineWidth', 0.5)
plot(Ax(ai), [bouts(selbout, 2), EEG.xmax.*2.95/4]./(24*60*60), [-3.9 -6], ':', 'Color', [0.8, 0.8, 0.8], 'LineWidth', 0.5)
% Time scale
plot(Ax(ai), [EEG.xmax-3600, EEG.xmax]./(24*60*60), [-4, -4], '-k', 'LineWidth', 2)
text(Ax(ai), mean([EEG.xmax-3600, EEG.xmax]./(24*60*60)), -4, '1h', ...
    'FontSize', 6, ...
    'VerticalAlignment', 'top', ...
    'HorizontalAlignment', 'center')
% Axes props
Ax(ai).YLim = [-3.75, 1.75];
Ax(ai).Box = 'off';
Ax(ai).XColor = 'w';
Ax(ai).TickLength = [0, 0];
Ax(ai).XTick = [];
Ax(ai).FontSize = 8;
Ax(ai).Clipping = 'off';

% -------------------------------------------------------------------------
% Number of bouts
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Create axes
ai = ai+1;
Ax(ai) = axes(...
    'NextPlot', 'add', ...
    'FontSize', 6, ...
    'Box', 'on', ...
    'OuterPosition', [8/12, 1-3/12, 1/12, 3/12]+[-0.5/12 0.02 0 -0.02], ...
    'YLim', [0, 18], ...
    'YTick', 0:1:18, ...
    'Color', [0.96 0.97 0.99]);
YData = [BOUTS.a.nbouts(strcmpi(BOUTS.a.ses, 'placebo')), BOUTS.a.nbouts(strcmpi(BOUTS.a.ses, 'etc120'))];
EData = std(YData);
YData = mean(YData);
b = bar(Ax(ai), 1, YData, 'LineWidth', 0.5);
b(1).FaceColor = standard_colors('green');
b(2).FaceColor = standard_colors('yellow');
errorbar(Ax(ai), [0.85 1.15], YData, EData, 'k', 'LineStyle', 'none', 'CapSize', 2);
Ax(ai).XTick = [0.85 1.15];
Ax(ai).YGrid = 'on';
Ax(ai).XTickLabel = {'PLB', 'ETC'};
Ax(ai).XTickLabelRotation = 90;
Ax(ai).YTickLabel(2:end-1) = {''};
Ax(ai).YLabel.String = 'N';
Ax(ai).YLabel.FontSize = 8;
Ax(ai).TickLength = [0, 0];

% -------------------------------------------------------------------------
% Bout duration
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Create axes
ai = ai+1;
Ax(ai) = axes(...
    'NextPlot', 'add', ...
    'FontSize', 6, ...
    'Box', 'on', ...
    'OuterPosition', [9/12, 1-3/12, 1/12, 3/12]+[-0.25/12 0.02 0.0049 -0.02], ...
    'YLim', [0, 1200], ...
    'YTick', 0:100:1200, ...
    'Color', [0.96 0.97 0.99]);
YData = {BOUTS.b.dur(strcmpi(BOUTS.b.ses, 'placebo')), BOUTS.b.dur(strcmpi(BOUTS.b.ses, 'etc120'))};
EData = cellfun(@(c) std(c), YData);
YData = cellfun(@(c) mean(c), YData);
b = bar(Ax(ai), 1, YData, 'LineWidth', 0.5);
b(1).FaceColor = standard_colors('green');
b(2).FaceColor = standard_colors('yellow');
errorbar(Ax(ai), [0.85 1.15], YData, EData, 'k', 'LineStyle', 'none', 'CapSize', 2);
Ax(ai).XTick = [0.85 1.15];
Ax(ai).YGrid = 'on';
Ax(ai).XTickLabel = {'PLB', 'ETC'};
Ax(ai).XTickLabelRotation = 90;
Ax(ai).YTickLabel(2:end-3) = {''};
Ax(ai).YTickLabel(end-1:end) = {''};
Ax(ai).YTickLabel(end-2) = {'10^{3}'};
Ax(ai).YLabel.String = 'duration (s)';
Ax(ai).YLabel.FontSize = 8;
Ax(ai).TickLength = [0, 0];

% -------------------------------------------------------------------------
% Bout proportion
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Create axes
ai = ai+1;
Ax(ai) = axes(...
    'NextPlot', 'add', ...
    'FontSize', 6, ...
    'Box', 'on', ...
    'XLim', [0, 100], ...
    'XTick', [0, 100], ...
    'YLim', [0, 1], ...
    'YTick', 0:0.1:1, ...
    'OuterPosition', [10/12, 1-3/12, 2/12, 3/12]+[0 0.0376 0 -0.0376], ...
    'Color', [0.96 0.97 0.99]);
histogram(Ax(ai), [BOUTS.b.onset(strcmpi(BOUTS.b.ses, 'placebo')).*100; 101], 100, ...
    'EdgeColor', standard_colors('green'), ...
    'LineWidth', 1, ...
    'Normalization', 'cdf', ...
    'DisplayStyle', 'stairs');
histogram(Ax(ai), [BOUTS.b.onset(strcmpi(BOUTS.b.ses, 'etc120')).*100; 101], 100, ...
    'EdgeColor', standard_colors('yellow'), ...
    'LineWidth', 1, ...
    'Normalization', 'cdf', ...
    'DisplayStyle', 'stairs');
Ax(ai).YGrid = 'on';
Ax(ai).YTickLabel(2:end-1) = {''};
Ax(ai).XLabel.String = 'onset (%)';
Ax(ai).XLabel.FontSize = 8;
Ax(ai).YLabel.String = 'cdf';
Ax(ai).YLabel.FontSize = 8;
Ax(ai).XLabel.Position(2) = -0.05;
Ax(ai).TickLength = [0, 0];

% -------------------------------------------------------------------------
% Placebo Spectra
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Create axes
ai = ai+1;
Ax(ai) = axes(...
    'NextPlot', 'add', ...
    'FontSize', 8, ...
    'XTick', 0:0.01:0.1, ...
    'XTickLabel', {'0', '', '0.02', '', '' , '', '', '', '', '', '0.1'}, ...
    'YLim', AmpYLim, ...
    'YTick', 1:10, ...
    'XGrid', 'on', ...
    'TickLength', [0 0], ...
    'Box', 'on', ...
    'OuterPosition', [6/12 1-6.5/12, 1.5/12, 4/12] + [-0.0353 0 0.25/12 -0.03], ...
    'Color', [0.96 0.97 0.99]);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Extact X and Y-data
idx_cond = strcmpi({ISF.session}, 'placebo');
XData = repmat([asrow(ISF(1).freqs), nan], 1, 178);
YData = cat(3, ISF(idx_cond).data);
YData = mean(YData, 3);
YData = [YData'; nan(1, 178)];
YData = ascolumn(YData(:));
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Plot as patch 
patch(Ax(ai), 'XData', XData, 'YData', YData, ...
    'EdgeColor', 'k', ...
    'EdgeAlpha', 0.05, ...
    'FaceColor', 'none');
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Set axes props
Ax(ai).FontSize = 6;
Ax(ai).Title.String = 'PLACEBO';
Ax(ai).Title.FontWeight = 'normal';
Ax(ai).Title.FontSize = 8;
Ax(ai).XLabel.String = 'frequency (Hz)';
Ax(ai).XLabel.FontSize = 8;
Ax(ai).XTickLabelRotation = 0;
Ax(ai).YLabel.String = 'amp. (a.u.)';
Ax(ai).YLabel.FontSize = 8;

% -------------------------------------------------------------------------
% ETC 120 Spectra
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Create axes
ai = ai+1;
Ax(ai) = axes(...
    'NextPlot', 'add', ...
    'FontSize', 8, ...
    'XTick', 0:0.01:0.1, ...
    'XTickLabel', {'0', '', '0.02', '', '' , '', '', '', '', '', '0.1'}, ...
    'YLim', AmpYLim, ...
    'YTick', 1:10, ...
    'XGrid', 'on', ...
    'TickLength', [0 0], ...
    'Box', 'on', ...
    'OuterPosition', [7.5/12 1-6.5/12, 1.5/12, 4/12] + [-0.25/12 0 0.25/12 -0.03], ...
    'Color', [0.96 0.97 0.99]);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Extact X and Y-data
idx_cond = strcmpi({ISF.session}, 'etc120');
XData = repmat([asrow(ISF(1).freqs), nan], 1, 178);
YData = cat(3, ISF(idx_cond).data);
YData = mean(YData, 3);
YData = [YData'; nan(1, 178)];
YData = ascolumn(YData(:));
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Plot as patch 
patch(Ax(ai), 'XData', XData, 'YData', YData, ...
    'EdgeColor', 'k', ...
    'EdgeAlpha', 0.05, ...
    'FaceColor', 'none');
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Set axes props
Ax(ai).FontSize = 6;
Ax(ai).Title.String = 'ETC120';
Ax(ai).Title.FontWeight = 'normal';
Ax(ai).Title.FontSize = 8;
Ax(ai).XLabel.String = 'frequency (Hz)';
Ax(ai).XLabel.FontSize = 8;
Ax(ai).XTickLabelRotation = 0;

% -------------------------------------------------------------------------
% PLACEBO X-Correlation (Sigma and HR)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Create axes
ai = ai+1;
Ax(ai) = axes(...
    'NextPlot', 'add', ...
    'FontSize', 8, ...
    'XTick', XC(1).xmin:10:XC(1).xmax, ...
    'YLim', [-0.15, 0.15], ...
    'YTick', -1:0.1:1, ...
    'XGrid', 'on', ...
    'TickLength', [0 0], ...
    'OuterPosition', [9/12 1-6.5/12, 1.5/12, 4/12] + [-0.25/12 0 0.25/12 -0.03], ...
    'Box', 'on', ...
    'Color', [0.96 0.97 0.99]);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Extact X and Y-data
idx_cond = strcmpi({XC.session}, 'placebo');
XData = repmat([asrow(XC(1).times), nan], 1, 178);
YData = cat(3, XC(idx_cond).data);
YData = mean(YData, 3);
[~, idx_peak] = max(YData');
mlag = mean(XC(1).times(idx_peak));
YData = [YData'; nan(1, 178)];
YData = ascolumn(YData(:));
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Plot as patch 
patch(Ax(ai), 'XData', XData, 'YData', YData, ...
    'EdgeColor', 'k', ...
    'EdgeAlpha', 0.05, ...
    'FaceColor', 'none');
plot(Ax(ai), [mlag, mlag], [-0.15 0.15], ':k')
text(Ax(ai), mlag, -0.1, sprintf(' %.1f s', mlag), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'baseline', 'FontSize', 6)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Set axes props
Ax(ai).FontSize = 6;
Ax(ai).Title.String = 'PLACEBO';
Ax(ai).Title.FontWeight = 'normal';
Ax(ai).Title.FontSize = 8;
Ax(ai).XLabel.String = 'time lag (s)';
Ax(ai).XLabel.FontSize = 8;
Ax(ai).XTickLabelRotation = 0;
Ax(ai).YLabel.String = 'r';
Ax(ai).YLabel.FontSize = 8;
Ax(ai).XTickLabel(2:end-(length(Ax(ai).XTickLabel)-1)/2-1) = {''};
Ax(ai).XTickLabel(end-(length(Ax(ai).XTickLabel)-1)/2+1:end-1) = {''};
tmp = Ax(ai).Position;
Ax(ai).YLabel.Position(1) = XC(1).xmin-25;
Ax(ai).Position = tmp;

% -------------------------------------------------------------------------
% ETC120 X-Correlation (Sigma and HR)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Create axes
ai = ai+1;
Ax(ai) = axes(...
    'NextPlot', 'add', ...
    'FontSize', 8, ...
    'XTick', XC(1).xmin:10:XC(1).xmax, ...
    'YLim', [-0.15, 0.15], ...
    'YTick', -1:0.1:1, ...
    'XGrid', 'on', ...
    'TickLength', [0 0], ...
    'OuterPosition', [10.5/12 1-6.5/12, 1.5/12, 4/12] + [-0.25/12 0 0.25/12 -0.03], ...
    'Box', 'on', ...
    'Color', [0.96 0.97 0.99]);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Extact X and Y-data
idx_cond = strcmpi({XC.session}, 'etc120');
XData = repmat([asrow(XC(1).times), nan], 1, 178);
YData = cat(3, XC(idx_cond).data);
YData = mean(YData, 3);
[~, idx_peak] = max(YData');
mlag = mean(XC(1).times(idx_peak));
YData = [YData'; nan(1, 178)];
YData = ascolumn(YData(:));
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Plot as patch 
patch(Ax(ai), 'XData', XData, 'YData', YData, ...
    'EdgeColor', 'k', ...
    'EdgeAlpha', 0.05, ...
    'FaceColor', 'none');
plot(Ax(ai), [mlag, mlag], [-0.15 0.15], ':k')
text(Ax(ai), mlag, -0.1, sprintf(' %.1f s', mlag), 'HorizontalAlignment', 'left', 'VerticalAlignment', 'baseline', 'FontSize', 6)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Set axes props
Ax(ai).FontSize = 6;
Ax(ai).Title.String = 'ETC120';
Ax(ai).Title.FontWeight = 'normal';
Ax(ai).Title.FontSize = 8;
Ax(ai).XLabel.String = 'time lag (s)';
Ax(ai).XLabel.FontSize = 8;
Ax(ai).XTickLabelRotation = 0;
Ax(ai).YLabel.String = 'r';
Ax(ai).YLabel.FontSize = 8;
Ax(ai).XTickLabel(2:end-(length(Ax(ai).XTickLabel)-1)/2-1) = {''};
Ax(ai).XTickLabel(end-(length(Ax(ai).XTickLabel)-1)/2+1:end-1) = {''};
tmp = Ax(ai).Position;
Ax(ai).YLabel.Position(1) = XC(1).xmin-25;
Ax(ai).Position = tmp;

% -------------------------------------------------------------------------
% ISO FREQUENCY TOPOPLOTS
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Create background
ai = ai+1;
Ax(ai) = axes('Position', [0/9 1-12/12 3/9 5/12]+[1/48 1/48 -1/24 -1/24], ...
    'Box', 'on', ...
    'XColor', [0.86 0.87 0.89], ...
    'YColor', [0.86 0.87 0.89], ...
    'XTick', [], ...
    'YTick', [], ...
    'Color', [0.96 0.97 0.99]);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Create axes
ai = ai+1;
Ax(ai) = axes('Position', [0/9 1-12/12 1/9 5/12] + [0.03 -0.06 -0.02 0]);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Extact Y-data
idx_cond = strcmpi({ISF.session}, 'placebo');
YData = arrayfun(@(s) s.features(1).data, ISF(idx_cond), 'UniformOutput', false);
YData = cat(2, YData{:});
YData = mean(YData, 2);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Plot as topoplot 
topoplot(YData, chanlocs, ...
    'colormap', batlow, ...
    'hlinewidth', 1, ...
    'maplimits', [0.017, 0.023], ...
    'contourvals', (YData > 0.019) + (YData > 0.021), ...
    'numcontour', length(unique((YData > 0.019) + (YData > 0.021)))-1, ...
    'conv', 'on', ...
    'whitebk', 'on');
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Ax(ai).Title.String = 'PLACEBO';
Ax(ai).Title.FontSize = 8;
Ax(ai).Title.FontWeight = 'normal';
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% SET COLORBAR
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
h = colorbar(Ax(ai));
h.Location = 'northoutside';
h.Position = [1/12 1-8/12 1/12 1/48] + [-0.02 -0.02 0.04 0];
h.Ticks = [0.017, 0.019, 0.021, 0.023];
h.TickLength = 0.07;
h.FontSize = 6;
h.Label.String = 'Sigma ISF frequency (Hz)';
h.Label.FontSize = 8;
h.Label.FontWeight = 'bold';
h.Label.Position(2) = 4.2;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Create axes
ai = ai+1;
Ax(ai) = axes('Position', [1/9 1-12/12 1/9 5/12] + [0.01 -0.06 -0.02 0]);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Extact Y-data
idx_cond = strcmpi({ISF.session}, 'etc120');
YData = arrayfun(@(s) s.features(1).data, ISF(idx_cond), 'UniformOutput', false);
YData = cat(2, YData{:});
YData = mean(YData, 2);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Plot as topoplot 
topoplot(YData, chanlocs, ...
    'colormap', batlow, ...
    'hlinewidth', 1, ...
    'maplimits', [0.017, 0.023], ...
    'contourvals', (YData > 0.019) + (YData > 0.021), ...
    'numcontour', length(unique((YData > 0.019) + (YData > 0.021)))-1, ...
    'conv', 'on', ...
    'whitebk', 'on');
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Ax(ai).Title.String = 'ETC120';
Ax(ai).Title.FontSize = 8;
Ax(ai).Title.FontWeight = 'normal';
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Create axes
ai = ai+1;
Ax(ai) = axes('Position', [2/9 1-12/12 1/9 5/12] + [-0.01 -0.06 -0.02 0]);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Extact Y-data and P-value data
YData = GLM.result(1).t.stat;
PData = GLM.result(1).t.p_fwe;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Plot as topoplot 
topoplot(YData, chanlocs, ...
    'colormap', batlow, ...
    'hlinewidth', 1, ...
    'maplimits', [-3.2905 3.2905], ...
    'emarker2', {find(PData < 0.05), '.', 'w', 8, 1}, ... 
    'contourvals', (YData > -1.9623) + (YData > 1.9623), ...
    'numcontour', length(unique((YData > -1.9623) + (YData > 1.9623)))-1, ...
    'conv', 'on', ...
    'whitebk', 'on');
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Ax(ai).Title.String = 't-stat';
Ax(ai).Title.FontSize = 8;
Ax(ai).Title.FontWeight = 'normal';
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% SET COLORBAR
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
h = colorbar(Ax(ai));
h.Location = 'northoutside';
h.Position = [3/12 1-8/12 1/24 1/48] + [-0.035 -0.02 0.04 0];
h.Ticks = [-3.2905, -1.9623, 1.9623, 3.2905];
h.TickLabels = {'' '-1.96', '1.96', ''};
h.TickLength = 0.1;
h.FontSize = 6;
h.Ruler.TickLabelRotation = 0;

% -------------------------------------------------------------------------
% ISO PEAK AMP TOPOPLOTS
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Create background
ai = ai+1;
Ax(ai) = axes('Position', [3/9 1-12/12 3/9 5/12]+[1/48 1/48 -1/24 -1/24], ...
    'Box', 'on', ...
    'XColor', [0.86 0.87 0.89], ...
    'YColor', [0.86 0.87 0.89], ...
    'XTick', [], ...
    'YTick', [], ...
    'Color', [0.96 0.97 0.99]);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Create axes
ai = ai+1;
Ax(ai) = axes('Position', [3/9 1-12/12 1/9 5/12] + [0.03 -0.06 -0.02 0]);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Extact Y-data
idx_cond = strcmpi({ISF.session}, 'placebo');
YData = arrayfun(@(s) s.features(2).data, ISF(idx_cond), 'UniformOutput', false);
YData = cat(2, YData{:});
YData = mean(YData, 2);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Plot as topoplot 
topoplot(YData, chanlocs, ...
    'colormap', batlow, ...
    'hlinewidth', 1, ...
    'maplimits', AmpYLim, ...
    'contourvals', (YData > 1) + (YData > 2) + (YData > 3) + (YData > 4), ...
    'numcontour', length(unique((YData > 1) + (YData > 2) + (YData > 3) + (YData > 4)))-1, ...
    'conv', 'on', ...
    'whitebk', 'on');
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Ax(ai).Title.String = 'PLACEBO';
Ax(ai).Title.FontSize = 8;
Ax(ai).Title.FontWeight = 'normal';
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% SET COLORBAR
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
h = colorbar(Ax(ai));
h.Location = 'northoutside';
h.Position = [5/12 1-8/12 1/12 1/48] + [-0.02 -0.02 0.04 0];
h.Ticks = 0:4;
h.TickLength = 0.07;
h.FontSize = 6;
h.Label.String = 'Sigma ISF Amplitude (a.u.)';
h.Label.FontSize = 8;
h.Label.FontWeight = 'bold';
h.Label.Position(2) = 4.2;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Create axes
ai = ai+1;
Ax(ai) = axes('Position', [4/9 1-12/12 1/9 5/12] + [0.01 -0.06 -0.02 0]);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Extact Y-data
idx_cond = strcmpi({ISF.session}, 'etc120');
YData = arrayfun(@(s) s.features(2).data, ISF(idx_cond), 'UniformOutput', false);
YData = cat(2, YData{:});
YData = mean(YData, 2);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Plot as topoplot 
topoplot(YData, chanlocs, ...
    'colormap', batlow, ...
    'hlinewidth', 1, ...
    'maplimits', AmpYLim, ...
    'contourvals', (YData > 1) + (YData > 2) + (YData > 3) + (YData > 4), ...
    'numcontour', length(unique((YData > 1) + (YData > 2) + (YData > 3) + (YData > 4)))-1, ...
    'conv', 'on', ...
    'whitebk', 'on');
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Ax(ai).Title.String = 'ETC120';
Ax(ai).Title.FontSize = 8;
Ax(ai).Title.FontWeight = 'normal';
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Create axes
ai = ai+1;
Ax(ai) = axes('Position', [5/9 1-12/12 1/9 5/12] + [-0.01 -0.06 -0.02 0]);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Extact Y-data and P-value data
YData = GLM.result(2).t.stat;
PData = GLM.result(2).t.p_fwe;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Plot as topoplot 
topoplot(YData, chanlocs, ...
    'colormap', batlow, ...
    'hlinewidth', 1, ...
    'maplimits', [-3.2905 3.2905], ...
    'emarker2', {find(PData < 0.05), '.', 'w', 8, 1}, ... 
    'contourvals', (YData > -1.9623) + (YData > 1.9623), ...
    'numcontour', length(unique((YData > -1.9623) + (YData > 1.9623)))-1, ...
    'conv', 'on', ...
    'whitebk', 'on');
Ax(ai).Colormap = roma;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Ax(ai).Title.String = 't-stat';
Ax(ai).Title.FontSize = 8;
Ax(ai).Title.FontWeight = 'normal';
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% SET COLORBAR
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
h = colorbar(Ax(ai));
h.Location = 'northoutside';
h.Position = [7/12 1-8/12 1/24 1/48] + [-0.035 -0.02 0.04 0];
h.Ticks = [-3.2905, -1.9623, 1.9623, 3.2905];
h.TickLabels = {'' '-1.96', '1.96', ''};
h.TickLength = 0.1;
h.FontSize = 6;
h.Ruler.TickLabelRotation = 0;

% -------------------------------------------------------------------------
% XCORR LAG TOPOPLOTS (PLACEBO)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Create background
ai = ai+1;
Ax(ai) = axes('Position', [6/9 1-12/12 3/9 5/12]+[1/48 1/48 -1/24 -1/24], ...
    'Box', 'on', ...
    'XColor', [0.86 0.87 0.89], ...
    'YColor', [0.86 0.87 0.89], ...
    'XTick', [], ...
    'YTick', [], ...
    'Color', [0.96 0.97 0.99]);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Create axes
ai = ai+1;
Ax(ai) = axes('Position', [6/9 1-12/12 1/9 5/12] + [0.03 -0.06 -0.02 0]);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Extact Y-data
idx_cond = strcmpi({XC.session}, 'placebo');
YData = arrayfun(@(s) ascolumn(s.features(1).data), XC(idx_cond), 'UniformOutput', false);
YData = cat(2, YData{:});
YData = mean(YData, 2);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Plot as topoplot 
topoplot(YData, chanlocs, ...   
    'colormap', batlow, ...
    'hlinewidth', 1, ...
    'maplimits', [1, 4], ...
    'contourvals', (YData > 1) + (YData > 2) + (YData > 3) + (YData > 4), ...
    'numcontour', length(unique((YData > 1) + (YData > 2) + (YData > 3) + (YData > 4)))-1, ...
    'conv', 'on', ...
    'whitebk', 'on');
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Ax(ai).Title.String = 'PLACEBO';
Ax(ai).Title.FontSize = 8;
Ax(ai).Title.FontWeight = 'normal';
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% SET COLORBAR
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
h = colorbar(Ax(ai));
h.Limits = [0.8 4.2];
h.Location = 'northoutside';
h.Position = [9/12 1-8/12 1/12 1/48] + [-0.02 -0.02 0.04 0];
h.Ticks = 0:10;
h.TickLength = 0.07;
h.FontSize = 6;
h.Label.String = 'Time lag (s)';
h.Label.FontSize = 8;
h.Label.FontWeight = 'bold';
h.Label.Position(1) = 1;
h.Label.Position(2) = 4.2;

% -------------------------------------------------------------------------
% XCORR LAG TOPOPLOTS (ETC120)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Create axes
ai = ai+1;
Ax(ai) = axes('Position', [7/9 1-12/12 1/9 5/12] + [0.01 -0.06 -0.02 0]);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Extact Y-data
idx_cond = strcmpi({XC.session}, 'etc120');
YData = arrayfun(@(s) ascolumn(s.features(1).data), XC(idx_cond), 'UniformOutput', false);
YData = cat(2, YData{:});
YData = mean(YData, 2);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Plot as topoplot 
topoplot(YData, chanlocs, ...
    'colormap', batlow, ...
    'hlinewidth', 1, ...
    'maplimits', [1, 4], ...
    'contourvals', (YData > 1) + (YData > 2) + (YData > 3) + (YData > 4), ...
    'numcontour', length(unique((YData > 1) + (YData > 2) + (YData > 3) + (YData > 4)))-1, ...
    'conv', 'on', ...
    'whitebk', 'on');
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Ax(ai).Title.String = 'ETC120';
Ax(ai).Title.FontSize = 8;
Ax(ai).Title.FontWeight = 'normal';
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Create axes
ai = ai+1;
Ax(ai) = axes('Position', [8/9 1-12/12 1/9 5/12] + [-0.01 -0.06 -0.02 0]);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Extact Y-data and P-value data
YData = GLMXC.result(1).t.stat;
PData = GLMXC.result(1).t.p_fwe;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Plot as topoplot 
topoplot(YData, chanlocs, ...
    'colormap', batlow, ...
    'hlinewidth', 1, ...
    'maplimits', [-3.2905 3.2905], ...
    'emarker2', {find(PData < 0.05), '.', 'w', 8, 1}, ... 
    'contourvals', (YData > -1.9623) + (YData > 1.9623), ...
    'numcontour', length(unique((YData > -1.9623) + (YData > 1.9623)))-1, ...
    'conv', 'on', ...
    'whitebk', 'on');
Ax(ai).Colormap = roma;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Ax(ai).Title.String = 't-stat';
Ax(ai).Title.FontSize = 8;
Ax(ai).Title.FontWeight = 'normal';
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% SET COLORBAR
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
h = colorbar(Ax(ai));
h.Location = 'northoutside';
h.Position = [11/12 1-8/12 1/24 1/48] + [-0.035 -0.02 0.04 0];
h.Ticks = [-3.2905, -1.9623, 1.9623, 3.2905];
h.TickLabels = {'' '-1.96', '1.96', ''};
h.TickLength = 0.1;
h.FontSize = 6;
h.Ruler.TickLabelRotation = 0;

% -------------------------------------------------------------------------
% Sigma power timeseries
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
add = [540 -120];
% Create axes
ai = ai+1;
CLim = [0 32];
Ax(ai) = axes(...
    'NextPlot', 'add', ...
    'FontSize', 10, ...
    'CLim', CLim, ...
    'Box', 'on', ...
    'OuterPosition', [0/12, 1-6/12, 6/12, 3/12] + [0 0.05 0 0], ...
    'Color', [0.96 0.97 0.99]);
idx_bouts = find(strcmpi({SIGMA.event.type}, 'boundary'));
idx_smp = round([SIGMA.event(idx_bouts(selbout-1)).latency, SIGMA.event(idx_bouts(selbout)).latency])+add.*SIGMA.srate;
idx_chan = [chanlocs.cluster] > 0;
XData = ascolumn(SIGMA.times(idx_smp(1):idx_smp(2))./(24*60*60));
XData = [XData; XData(end)+1; -1];
YData = ascolumn(mean(SIGMA.data(idx_chan, idx_smp(1):idx_smp(2)), 1));
YData = [YData; -6; -6];
patch(Ax(ai), 'XData', XData, 'YData', YData, 'CData', YData, ...
    'EdgeColor', standard_colors('blue'), ...
    'FaceColor', 'interp')
CMap = standard_colors('blue');
Ax(ai).Colormap = [linspace(1, CMap(1), 256); linspace(1, CMap(2), 256); linspace(1, CMap(3), 256)]';
% YLabel
text(Ax(ai), SIGMA.times(idx_smp(1))./(24*60*60), mean(CLim), 0, '\sigma power ', ...
    'FontSize', 8, ...
    'VerticalAlignment', 'top', ...
    'HorizontalAlignment', 'right')
% Axes props
Ax(ai).XLim = SIGMA.times(idx_smp)./(24*60*60);
Ax(ai).YLim = CLim;
Ax(ai).XTick = [];
Ax(ai).YTick = [];
Ax(ai).Visible = 'off';

% -------------------------------------------------------------------------
% HR timeseries
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Create axes
ai = ai+1;
CLim = [60 120];
Ax(ai) = axes(...
    'NextPlot', 'add', ...
    'FontSize', 10, ...
    'Box', 'on', ...
    'CLim', CLim, ...
    'OuterPosition', [0/12, 1-7/12, 6/12, 2/12] + [0 0.08 0 0], ...
    'Color', [0.96 0.97 0.99]);
idx_smp = round((bouts(selbout, :)+add).*HR.srate);
XData = ascolumn(HR.times(idx_smp(1):idx_smp(2))./(24*60*60));
XData = [XData; XData(end)+1; -1];
YData = ascolumn(mean(HR.data(1, idx_smp(1):idx_smp(2)), 1));
YData = [YData; 0; 0];
patch(Ax(ai), 'XData', XData, 'YData', YData, 'CData', YData, ...
    'EdgeColor', standard_colors('turquoise'), ...
    'FaceColor', 'interp')
CMap = standard_colors('turquoise');
Ax(ai).Colormap = [linspace(1, CMap(1), 256); linspace(1, CMap(2), 256); linspace(1, CMap(3), 256)]';
% Time scale
plot(Ax(ai), [max(HR.times(idx_smp))-50, max(HR.times(idx_smp))]./(24*60*60), [CLim(1) CLim(1)], '-k', 'LineWidth', 2)
text(Ax(ai), mean([max(HR.times(idx_smp))-50, max(HR.times(idx_smp))]./(24*60*60)), CLim(1), '50s', ...
    'FontSize', 6, ...
    'VerticalAlignment', 'top', ...
    'HorizontalAlignment', 'center')
% YLabel
text(Ax(ai), HR.times(idx_smp(1))./(24*60*60), mean(CLim), 0, 'HR ', ...
    'FontSize', 8, ...
    'VerticalAlignment', 'top', ...
    'HorizontalAlignment', 'right')
% Axes props
Ax(ai).XLim = HR.times(idx_smp)./(24*60*60);
Ax(ai).YLim = CLim;
Ax(ai).XTick = [];
Ax(ai).YTick = [];
Ax(ai).Visible = 'off';

Ax(ai-10).Colormap = roma;
Ax(ai-6).Colormap = roma;

% -------------------------------------------------------------------------
% Box labels
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Create axes
ai = ai+1;
Ax(ai) = axes('Position', [0 0 1 1], 'XLim', [0 12], 'YLim', [0 12]);
Ax(ai).Visible = 'off';
text(Ax(ai), 0, 12, 'A', 'FontSize', 8, 'FontWeight', 'bold', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')
text(Ax(ai), 7.5, 12, 'B', 'FontSize', 8, 'FontWeight', 'bold', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')
text(Ax(ai), 0, 9, 'C', 'FontSize', 8, 'FontWeight', 'bold', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')
text(Ax(ai), 5.6, 9.17, 'D', 'FontSize', 8, 'FontWeight', 'bold', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')
text(Ax(ai), 8.92, 9.17, 'E', 'FontSize', 8, 'FontWeight', 'bold', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')
text(Ax(ai), 0, 5.4, 'F', 'FontSize', 8, 'FontWeight', 'bold', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')
text(Ax(ai), 4, 5.4, 'G', 'FontSize', 8, 'FontWeight', 'bold', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')
text(Ax(ai), 8, 5.4, 'H', 'FontSize', 8, 'FontWeight', 'bold', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top')

disp('Saving figure...')
exportgraphics(Fig, sprintf('figures/fig_%sisffeatures.png', type), 'Resolution', 300)
disp('Done saving')
