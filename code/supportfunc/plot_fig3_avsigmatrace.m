function [Ax, l] = plot_fig3_avsigmatrace(Fig, H, pcfg)

% Create axes
clear Ax
Ax = axes(Fig, 'NextPlot', 'add');
Ax.Color = [0.96 0.97 0.99];
Ax.Box = 'on';
Ax.FontSize = 8;
Ax.TickLength = [0, 0];
Ax.XLim = [-95, 65];
Ax.XTick = -90:30:60;
Ax.YLim = [-1.5 0.75];
Ax.YTick = [0, 0.75];
Ax.XGrid = 'on';
Ax.YGrid = 'on';

Ax.XLabel.String = 'time to REM (s)';
Ax.XLabel.FontSize = 10;
Ax.XTickLabelRotation = 0;

Ax.YLabel.String = '\sigma {\fontsize{8}(dB)}';
Ax.YLabel.Interpreter = 'tex';
Ax.YLabel.FontSize = 12;
Ax.YLabel.Color = 'k';

clear idx YData EData
idx.pbo = find(strcmpi({H.cond}, 'placebo'));
idx.etc = find(strcmpi({H.cond}, 'etc120'));
YData.pbo = arrayfun(@(r) mean(smoothdata(H(r).rawsigma(2:3601, :), 'movmean', 100), 2), idx.pbo, 'UniformOutput', false)';
YData.etc = arrayfun(@(r) mean(smoothdata(H(r).rawsigma(2:3601, :), 'movmean', 100), 2), idx.etc, 'UniformOutput', false)';
YData.pbo = double(cat(2, YData.pbo{:}));
YData.etc = double(cat(2, YData.etc{:}));
XData = -300:0.1:59.9;

EData.pbo = tinv(0.975, size(YData.pbo, 2)-1).*(std(YData.pbo, [], 2)./sqrt(size(YData.pbo, 2)));
EData.etc = tinv(0.975, size(YData.etc, 2)-1).*(std(YData.etc, [], 2)./sqrt(size(YData.etc, 2)));

errorpatch(Ax, XData, mean(YData.pbo, 2), EData.pbo, ...
    'FaceColor', css_standard_colors('pbo'), ...
    'FaceAlpha', 0.3);

errorpatch(Ax, XData, mean(YData.etc, 2), EData.etc, ...
    'FaceColor', css_standard_colors('etc'), ...
    'FaceAlpha', 0.3);


clear p
p(1) = plot(Ax, XData, mean(YData.pbo, 2), '-', ...
    'LineWidth', 1, ...
    'Color', css_standard_colors('pbo'));
p(2) = plot(Ax, XData, mean(YData.etc, 2), '-', ...
    'LineWidth', 1, ...
    'Color', css_standard_colors('etc'));

plot(Ax, [0, 0], Ax.YLim, ':k', 'LineWidth', 1)

% Plot text
text(Ax, 6, Ax.YLim(2), ' N2 → REM ', ...
    'FontSize', 8, ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'top', ...
    'Color', 'k')

l = legend(p, {'placebo', 'THC/CBD'}, ...
    'Box', 'off', ...
    'FontSize', 8, ...
    'Location', 'northoutside');
l.ItemTokenSize = [10,10];

Ax.YLabel.Position(1) = -110;
end