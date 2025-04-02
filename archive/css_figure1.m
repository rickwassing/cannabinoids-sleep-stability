function css_figure1()

% Load group-level results, colormap, and chanlocs
load('group-level/unpaired-test_clustmass_perm/glm.mat')
load('colormap_roma.mat')
chanlocs = template_to_chanlocs(which('GSN-HydroCel-257.sfp'));
[~, system] = ishdeeg({chanlocs.labels});
incl = hdeeg_scalpchannels(system);
chanlocs = chanlocs(ismember({chanlocs.labels}, incl));

% Load data from both groups
clear CNT MCI
Files = dir('derivatives/EEG-output-fstlvl/sub-*/ses-*/sub-*_ismfit.mat');
PHEN = readtable('participants.csv');
cnt = struct();
cnt.cnt = 0;
cnt.mci = 0;
for i = 1:length(Files)
    kv = filename2struct(Files(i).name);
    idx_phen = strcmpi(PHEN.participant_id, kv.sub);
    switch PHEN.group{idx_phen}
        case 'control'
            cnt.cnt = cnt.cnt+1;
            CNT(cnt.cnt) = LoadDataset(fullfile(Files(i).folder, Files(i).name), 'matrix');
        case 'mci'
            cnt.mci = cnt.mci+1;
            MCI(cnt.mci) = LoadDataset(fullfile(Files(i).folder, Files(i).name), 'matrix');
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
% Control spectra
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Create axes
Ax(1) = axes(...
    'NextPlot', 'add', ...
    'FontSize', 10, ...
    'XTick', [0, 0.02, 0.04, 0.06, 0.1], ...
    'YLim', [0, 4], ...
    'YTick', 1:10, ...
    'XGrid', 'on', ...
    'TickLength', [0 0], ...
    'Box', 'on', ...
    'OuterPosition', [0, 0.78, 0.5, 0.23], ...
    'Color', [0.96 0.97 0.99]);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Extact X and Y-data
XData = repmat([CNT(1).freqs, nan], 1, 178);
YData = mean(cat(3, CNT.data), 3);
YData = [YData'; nan(1, 178)];
YData = YData(:);
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
% Ax(1).XLabel.Position(2) = -0.1;
Ax(1).YLabel.String = 'amplitude (a.u.)';
Ax(1).YLabel.FontSize = 8;
Ax(1).XTickLabelRotation = 0;

% -------------------------------------------------------------------------
% MCI spectra
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Create axes
Ax(2) = axes(...
    'NextPlot', 'add', ...
    'FontSize', 10, ...
    'XTick', [0, 0.02, 0.04, 0.06, 0.1], ...
    'YLim', [0, 4], ...
    'YTick', 1:10, ...
    'XGrid', 'on', ...
    'TickLength', [0 0], ...
    'Box', 'on', ...
    'OuterPosition', [0.5, 0.78, 0.5, 0.23], ...
    'Color', [0.96 0.97 0.99]);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Extact X and Y-data
XData = repmat([MCI(1).freqs, nan], 1, 178);
YData = mean(cat(3, MCI.data), 3);
YData = [YData'; nan(1, 178)];
YData = YData(:);
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
% Ax(2).XLabel.Position(2) = -0.1;
Ax(2).YLabel.String = 'amplitude (a.u.)';
Ax(2).YLabel.FontSize = 8;
Ax(2).XTickLabelRotation = 0;

% -------------------------------------------------------------------------
% CONTROL ISO FREQUENCY TOPOPLOTS
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Create axes
Ax(3) = axes('Position', [0, 0.3, 1/3, 0.6]);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Extact Y-data
YData = arrayfun(@(s) s.features(1).data, CNT, 'UniformOutput', false);
YData = cat(2, YData{:});
YData = mean(YData, 2);
tmpa = YData;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Plot as topoplot 
topoplot(YData, chanlocs, ...
    'contourvals', (YData > 0.019) + (YData > 0.021), ...
    'numcontour', length(unique((YData > 0.019) + (YData > 0.021)))-1, ...
    'conv', 'on', ...
    'whitebk', 'on');

% -------------------------------------------------------------------------
% MCI ISO FREQUENCY TOPOPLOTS
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Create axes
Ax(4) = axes('Position', [1/3, 0.3, 1/3, 0.6]);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Extact Y-data
YData = arrayfun(@(s) s.features(1).data, MCI, 'UniformOutput', false);
YData = cat(2, YData{:});
YData = mean(YData, 2);
tmpb = YData;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Plot as topoplot 
topoplot(YData, chanlocs, ...
    'contourvals', (YData > 0.019) + (YData > 0.021), ...
    'numcontour', length(unique((YData > 0.019) + (YData > 0.021)))-1, ...
    'conv', 'on', ...
    'whitebk', 'on');

% -------------------------------------------------------------------------
% CONTROL ISO AMPLITUDE TOPOPLOTS
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Create axes
Ax(5) = axes('Position', [0, -0.1, 1/3, 0.6]);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Extact Y-data
YData = arrayfun(@(s) s.features(2).data, CNT, 'UniformOutput', false);
YData = cat(2, YData{:});
YData = (mean(YData, 2));
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Plot as topoplot 
topoplot(YData, chanlocs, ...
    'contourvals', (YData > 2) + (YData > 3), ...
    'numcontour', length(unique((YData > 2) + (YData > 3)))-1, ...
    'conv', 'on', ...
    'whitebk', 'on');

% -------------------------------------------------------------------------
% MCI ISO AMPLITUDE TOPOPLOTS
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Create axes
Ax(6) = axes('Position', [1/3, -0.1, 1/3, 0.6]);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Extact Y-data
YData = arrayfun(@(s) s.features(2).data, MCI, 'UniformOutput', false);
YData = cat(2, YData{:});
YData = (mean(YData, 2));
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Plot as topoplot 
topoplot(YData, chanlocs, ...
    'contourvals', (YData > 2) + (YData > 3), ...
    'numcontour', length(unique((YData > 2) + (YData > 3)))-1, ...
    'conv', 'on', ...
    'whitebk', 'on');

% -------------------------------------------------------------------------
% CONTROL T-STAT TOPOPLOTS
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
% MCI T-STAT TOPOPLOTS
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
Ax(1).Title.String = 'CONTROLS';
Ax(2).Title.String = 'MCI';
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Ax(3).Colormap = roma;
Ax(3).CLim = [0.017, 0.023];
Ax(3).Title.String = 'CONTROLS';
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Ax(4).Colormap = roma;
Ax(4).CLim = [0.017, 0.023];
Ax(4).Title.String = 'MCI';
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Ax(5).Colormap = roma;
Ax(5).CLim = [1, 4];
Ax(5).Title.String = 'CONTROLS';
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Ax(6).Colormap = roma;
Ax(6).CLim = [1, 4];
Ax(6).Title.String = 'MCI';
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
h.Ticks = [1, 2, 3, 4];
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
exportgraphics(Fig, 'figures/fig_1.png', 'Resolution', 1200)
