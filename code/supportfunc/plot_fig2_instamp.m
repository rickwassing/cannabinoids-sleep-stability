function [Ax] = plot_fig2_instamp(Fig, T, pcfg)

% Create axes
Ax = axes(Fig, 'NextPlot', 'add');
Ax.Box = 'on';
Ax.FontSize = 8;
Ax.TickLength = [0 0];
Ax.YGrid = 'on';
Ax.Color = [0.96 0.97 0.99];

clear YData EData

cond = 'pbo';
fld = 'aw';
[~, YData.(cond).(fld)] = withinSubMean(T(pcfg.idx.(cond).(fld), :), 'amp_d0', 1);
YData.(cond).(fld) = YData.(cond).(fld)';
EData.(cond).(fld) = tinv(0.975, size(YData.(cond).(fld), 2)-1).*(std(YData.(cond).(fld), [], 2)./sqrt(size(YData.(cond).(fld), 2)));

cond = 'pbo';
fld = 'cs';
[~, YData.(cond).(fld)] = withinSubMean(T(pcfg.idx.(cond).(fld), :), 'amp_d0', 1);
YData.(cond).(fld) = YData.(cond).(fld)';
EData.(cond).(fld) = tinv(0.975, size(YData.(cond).(fld), 2)-1).*(std(YData.(cond).(fld), [], 2)./sqrt(size(YData.(cond).(fld), 2)));

cond = 'etc';
fld = 'aw';
[~, YData.(cond).(fld)] = withinSubMean(T(pcfg.idx.(cond).(fld), :), 'amp_d0', 1);
YData.(cond).(fld) = YData.(cond).(fld)';
EData.(cond).(fld) = tinv(0.975, size(YData.(cond).(fld), 2)-1).*(std(YData.(cond).(fld), [], 2)./sqrt(size(YData.(cond).(fld), 2)));

cond = 'etc';
fld = 'cs';
[~, YData.(cond).(fld)] = withinSubMean(T(pcfg.idx.(cond).(fld), :), 'amp_d0', 1);
YData.(cond).(fld) = YData.(cond).(fld)';
EData.(cond).(fld) = tinv(0.975, size(YData.(cond).(fld), 2)-1).*(std(YData.(cond).(fld), [], 2)./sqrt(size(YData.(cond).(fld), 2)));

h = errorbar(Ax, [1, 3], ...
    [mean(YData.pbo.aw), mean(YData.pbo.cs)], ...
    [EData.pbo.aw, EData.pbo.cs], 'o', ...
    'LineStyle', 'none', ...
    'Color', 'k', ...
    'MarkerSize', 3, ...
    'MarkerFaceColor', css_standard_colors('pbo'), ...
    'MarkerEdgeColor', css_standard_colors('pbo'));

h = errorbar(Ax, [2, 4], ...
    [mean(YData.etc.aw), mean(YData.etc.cs)], ...
    [EData.etc.aw, EData.etc.cs], 'o', ...
    'LineStyle', 'none', ...
    'Color', 'k', ...
    'MarkerSize', 3, ...
    'MarkerFaceColor', css_standard_colors('etc'), ...
    'MarkerEdgeColor', css_standard_colors('etc'));

Ax.XLim = [0.33, 4.67];
Ax.XTick = [1.5, 3.5];
Ax.XTickLabel = {'AW', 'CS'};
Ax.YLim = [0, 0.2];
Ax.YTick = [0, 0.2];
Ax.XLabel.String = '\sigma ISF Amp.';
Ax.XLabel.FontSize = 10;
Ax.YLabel.String = '\sigma {\fontsize{8}(dB^2)}';
Ax.YLabel.FontSize = 12;
Ax.YLabel.FontWeight = 'normal';
Ax.YLabel.Position(1) = Ax.XLim(1);

end