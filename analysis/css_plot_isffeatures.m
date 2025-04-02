function css_plot_isffeatures(type, varargin)

% Load group-level results, colormap, and chanlocs
load(sprintf('group-level/a1c_pairttest_cmass/%ssigma/glm.mat', type))
load('colormap_roma.mat')
load('colormap_batlow.mat')
chanlocs = template_to_chanlocs(which('GSN-HydroCel-257.sfp'));
[~, system] = ishdeeg({chanlocs.labels});
incl = hdeeg_scalpchannels(system);
chanlocs = chanlocs(ismember({chanlocs.labels}, incl));
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
% -------------------------------------------------------------------------
% Create new figure
close all
clear Ax
Fig = figure;
Fig.Units = 'centimeters';
Fig.Position = [5, 15, 8, 12];
% -------------------------------------------------------------------------
% Placebo Spectra
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Create axes
Ax(1) = axes(...
    'NextPlot', 'add', ...
    'FontSize', 10, ...
    'XTick', [0, 0.02, 0.04, 0.06, 0.1], ...
    'YLim', AmpYLim, ...
    'YTick', 1:10, ...
    'XGrid', 'on', ...
    'TickLength', [0 0], ...
    'Box', 'on', ...
    'OuterPosition', [0, 0.78, 0.5, 0.23], ...
    'Color', [0.96 0.97 0.99]);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Extact X and Y-data
idx_cond = strcmpi({ISF.session}, 'placebo');
XData = repmat([asrow(ISF(1).freqs), nan], 1, 178);
YData = cat(3, ISF(idx_cond).data);
YData = mean(YData(idx_chan, :, :), 3);
YData = [YData'; nan(1, length(idx_chan))];
YData = ascolumn(YData(:));
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Plot as patch 
patch(Ax(1), 'XData', XData, 'YData', YData, ...
    'EdgeColor', 'k', ...
    'EdgeAlpha', 0.05, ...
    'FaceColor', 'none');
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Set axes props
Ax(1).FontSize = 6;
Ax(1).Title.FontSize = 8;
Ax(1).XLabel.String = 'frequency (Hz)';
Ax(1).XLabel.FontSize = 8;
Ax(1).YLabel.String = 'amplitude (a.u.)';
Ax(1).YLabel.FontSize = 8;
Ax(1).XTickLabelRotation = 0;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if nargin > 1
    tmpax = axes();
    topoplot(rand(length(chanlocs), 1), chanlocs, ...
        'style', 'blank', ...
        'electrodes', 'on', ...
        'plotchans', idx_chan, ... 
        'whitebk', 'on')
    tmpax.Position = [0.28, 0.87, 0.20, 0.12];
    tmpax.XLim = [-1.1, 1.1];
    tmpax.YLim = [-1.1, 1.1];
    tmpax.Children(1).MarkerSize = 3;
    tmpax.Children(2).LineWidth = 0.5;
    tmpax.Children(3).LineWidth = 0.5;
    tmpax.Children(4).LineWidth = 0.5;
    tmpax.Children(5).LineWidth = 0.5;
end
% -------------------------------------------------------------------------
% ETC spectra
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Create axes
Ax(2) = axes(...
    'NextPlot', 'add', ...
    'FontSize', 10, ...
    'XTick', [0, 0.02, 0.04, 0.06, 0.1], ...
    'YLim', AmpYLim, ...
    'YTick', 1:10, ...
    'XGrid', 'on', ...
    'TickLength', [0 0], ...
    'Box', 'on', ...
    'OuterPosition', [0.5, 0.78, 0.5, 0.23], ...
    'Color', [0.96 0.97 0.99]);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Extact X and Y-data
idx_cond = strcmpi({ISF.session}, 'etc120');
XData = repmat([asrow(ISF(1).freqs), nan], 1, 178);
YData = cat(3, ISF(idx_cond).data);
YData = mean(YData(idx_chan, :, :), 3);
YData = [YData'; nan(1, length(idx_chan))];
YData = ascolumn(YData(:));
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Plot as patch 
patch(Ax(2), 'XData', XData, 'YData', YData, ...
    'EdgeColor', 'k', ...
    'EdgeAlpha', 0.05, ...
    'FaceColor', 'none');
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Set axes props
Ax(2).FontSize = 6;
Ax(2).Title.FontSize = 8;
Ax(2).XLabel.String = 'frequency (Hz)';
Ax(2).XLabel.FontSize = 8;
Ax(2).YLabel.String = 'amplitude (a.u.)';
Ax(2).YLabel.FontSize = 8;
Ax(2).XTickLabelRotation = 0;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if nargin > 1
    tmpax = axes();
    topoplot(rand(length(chanlocs), 1), chanlocs, ...
        'style', 'blank', ...
        'electrodes', 'on', ...
        'plotchans', idx_chan, ... 
        'whitebk', 'on')
    tmpax.Position = [0.78, 0.87, 0.20, 0.12];
    tmpax.XLim = [-1.1, 1.1];
    tmpax.YLim = [-1.1, 1.1];
    tmpax.Children(1).MarkerSize = 3;
    tmpax.Children(2).LineWidth = 0.5;
    tmpax.Children(3).LineWidth = 0.5;
    tmpax.Children(4).LineWidth = 0.5;
    tmpax.Children(5).LineWidth = 0.5;
end
% -------------------------------------------------------------------------
% PLACEBO ISO FREQUENCY TOPOPLOTS
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Create axes
Ax(3) = axes('Position', [0, 0.3, 1/3, 0.6]);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Extact Y-data
idx_cond = strcmpi({ISF.session}, 'placebo');
YData = arrayfun(@(s) s.features(1).data, ISF(idx_cond), 'UniformOutput', false);
YData = cat(2, YData{:});
YData = mean(YData, 2);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Plot as topoplot 
topoplot(YData, chanlocs, ...
    'contourvals', (YData > 0.019) + (YData > 0.021), ...
    'numcontour', length(unique((YData > 0.019) + (YData > 0.021)))-1, ...
    'conv', 'on', ...
    'whitebk', 'on');

% -------------------------------------------------------------------------
% ETC ISO FREQUENCY TOPOPLOTS
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Create axes
Ax(4) = axes('Position', [1/3, 0.3, 1/3, 0.6]);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Extact Y-data
idx_cond = strcmpi({ISF.session}, 'etc120');
YData = arrayfun(@(s) s.features(1).data, ISF(idx_cond), 'UniformOutput', false);
YData = cat(2, YData{:});
YData = mean(YData, 2);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Plot as topoplot 
topoplot(YData, chanlocs, ...
    'contourvals', (YData > 0.019) + (YData > 0.021), ...
    'numcontour', length(unique((YData > 0.019) + (YData > 0.021)))-1, ...
    'conv', 'on', ...
    'whitebk', 'on');

% -------------------------------------------------------------------------
% PLACEBO ISO AMPLITUDE TOPOPLOTS
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Create axes
Ax(5) = axes('Position', [0, -0.1, 1/3, 0.6]);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Extact Y-data
idx_cond = strcmpi({ISF.session}, 'placebo');
YData = arrayfun(@(s) s.features(2).data, ISF(idx_cond), 'UniformOutput', false);
YData = cat(2, YData{:});
YData = mean(YData, 2);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Plot as topoplot 
topoplot(YData, chanlocs, ...
    'contourvals', (YData > 1) + (YData > 2) + (YData > 3) + (YData > 4), ...
    'numcontour', length(unique((YData > 1) + (YData > 2) + (YData > 3) + (YData > 4)))-1, ...
    'conv', 'on', ...
    'whitebk', 'on');

% -------------------------------------------------------------------------
% ETC ISO AMPLITUDE TOPOPLOTS
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Create axes
Ax(6) = axes('Position', [1/3, -0.1, 1/3, 0.6]);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Extact Y-data
idx_cond = strcmpi({ISF.session}, 'etc120');
YData = arrayfun(@(s) s.features(2).data, ISF(idx_cond), 'UniformOutput', false);
YData = cat(2, YData{:});
YData = mean(YData, 2);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Plot as topoplot 
topoplot(YData, chanlocs, ...
    'contourvals', (YData > 1) + (YData > 2) + (YData > 3) + (YData > 4), ...
    'numcontour', length(unique((YData > 1) + (YData > 2) + (YData > 3) + (YData > 4)))-1, ...
    'conv', 'on', ...
    'whitebk', 'on');

% -------------------------------------------------------------------------
% ISF FREQUENCY T-STAT TOPOPLOTS
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Create axes
Ax(7) = axes('Position', [2/3, 0.3, 1/3, 0.6]);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Extact Y-data and P-value data
YData = GLM.result(1).t.stat;
PData = GLM.result(1).t.p_fwe;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Plot as topoplot 
topoplot(YData, chanlocs, ...
    'emarker2', {find(PData < 0.05), '.', 'w', 8, 1}, ... 
    'contourvals', (YData > -1.9623) + (YData > 1.9623), ...
    'numcontour', length(unique((YData > -1.9623) + (YData > 1.9623)))-1, ...
    'conv', 'on', ...
    'whitebk', 'on');

% -------------------------------------------------------------------------
% AMPLITUDE T-STAT TOPOPLOTS
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Create axes
Ax(8) = axes('Position', [2/3, -0.1, 1/3, 0.6]);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Extact Y-data and P-value data
YData = GLM.result(2).t.stat;
PData = GLM.result(2).t.p_fwe;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Plot as topoplot 
topoplot(YData, chanlocs, ...
    'emarker2', {find(PData < 0.05), '.', 'w', 8, 1}, ... 
    'contourvals', (YData > -1.9623) + (YData > 1.9623), ...
    'numcontour', length(unique((YData > -1.9623) + (YData > 1.9623)))-1, ...
    'conv', 'on', ...
    'whitebk', 'on');

% -------------------------------------------------------------------------
% SET MORE AXES PROPS
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Ax(1).Title.String = 'PLACEBO';
Ax(2).Title.String = 'ETC120';
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Ax(3).Colormap = batlow;
Ax(3).CLim = [0.017, 0.023];
Ax(3).Title.String = 'PLACEBO';
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Ax(4).Colormap = batlow;
Ax(4).CLim = [0.017, 0.023];
Ax(4).Title.String = 'ETC120';
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Ax(5).Colormap = batlow;
Ax(5).CLim = AmpYLim;
Ax(5).Title.String = 'PLACEBO';
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Ax(6).Colormap = batlow;
Ax(6).CLim = AmpYLim;
Ax(6).Title.String = 'ETC120';
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Ax(7).Colormap = roma;
Ax(7).CLim = [-3.2905 3.2905];
Ax(7).Title.String = ' ';
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Ax(8).Colormap = roma;
Ax(8).CLim = [-3.2905 3.2905];
Ax(8).Title.String = ' ';

% -------------------------------------------------------------------------
% SET COLORMAPS
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
h = colorbar(Ax(3));
h.Location = 'southoutside';
h.Position(1) = 1/6;
h.Position(2) = 43/100;
h.Position(3) = 1/3;
h.Ticks = [0.017, 0.019, 0.021, 0.023];
h.TickLength = 0.1;
h.FontSize = 6;
h.Label.String = 'Sigma ISO frequency (Hz)';
h.Label.FontSize = 8;
h.Label.FontWeight = 'bold';
h.Label.Position(2) = 2.5;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
h = colorbar(Ax(5));
h.Location = 'southoutside';
h.Position(1) = 1/6;
h.Position(2) = 3/100;
h.Position(3) = 1/3;
h.Ticks = 0:1:max(AmpYLim);
h.TickLength = 0.1;
h.FontSize = 6;
h.Label.String = 'Sigma ISO amplitude (a.u.)';
h.Label.FontSize = 8;
h.Label.FontWeight = 'bold';
h.Label.Position(2) = 2.5;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
h = colorbar(Ax(7));
h.Location = 'eastoutside';
h.Position(1) = 81.26042/100;
h.Position(2) = 34/100;
h.Position(3) = 1/24;
h.Ticks = [-3.2905, -1.9623, 1.9623, 3.2905];
h.TickLength = 0.225;
h.TickLabels = {'-3.29' '-1.96', '1.96', '3.29'};
h.FontSize = 6;
h.Label.String = 'T-stat';
h.Label.FontSize = 8;
h.Label.Position(1) = 3.4;

% -------------------------------------------------------------------------
% SET MORE AXES PROPS
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% These had to be set after the colorbars for some weird reason
Ax(3).Title.FontSize = 8;
Ax(4).Title.FontSize = 8;
Ax(5).Title.FontSize = 8;
Ax(6).Title.FontSize = 8;
Ax(7).Title.FontSize = 8;
Ax(8).Title.FontSize = 8;

% -------------------------------------------------------------------------
% SAVE FIGURE
disp('Saving figure...')
exportgraphics(Fig, sprintf('figures/fig_%sisffeatures.png', type), 'Resolution', 1200)
disp('Done saving')