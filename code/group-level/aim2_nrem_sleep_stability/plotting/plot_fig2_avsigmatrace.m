function [Ax] = plot_fig2_avsigmatrace(Fig, T, perms, fld, pcfg)

% Create axes
Ax = axes(Fig, 'NextPlot', 'add');
Ax.XLim = [-60, 5];
Ax.YLim = pcfg.ylim_pwr;
Ax.XTick = -60:10:10;
Ax.XGrid = 'on';
Ax.Color = [0.96 0.97 0.99];
Ax.Box = 'on';
Ax.YTick = [];
Ax.FontSize = 8;
Ax.XLabel.String = 'time to arousal (s)';
Ax.XLabel.FontSize = 10;
Ax.TickLength = [0, 0];
if strcmpi(fld, 'aw')
    Ax.YLabel.String = '\sigma {\fontsize{8}(dB)}';
    Ax.YLabel.Interpreter = 'tex';
    Ax.YLabel.FontSize = 12;
    Ax.YLabel.Color = 'k';
else
    Ax.YLabel.String = ' ';
    Ax.YLabel.FontSize = 12;
end

clear p

% Plot zero line and arousal marker
plot(Ax, Ax.XLim, [0 0], '-k')
plot(Ax, [0, 0], pcfg.ylim_pwr, ':k')

% Plot averaged sigma power
XData = (-59.9:0.1:0)-(abs(pcfg.const.crop_aro(2))/10-30);

% Plot significant bins
YData = double(perms(1).y.(fld).pval < 0.05);
YData(~YData) = nan;
plot(Ax, XData, YData*(Ax.YLim(1).*0.5), '-', 'LineWidth', 3, 'Color', [0.52 0.51 0.55]);

if isfield(pcfg.sigma_perm.(fld), 'perm_p')
    idx_pcor = find(pcfg.sigma_perm.(fld).perm_p < 0.05);
    for i = 1:length(idx_pcor)
        plot(Ax, XData(asrow(pcfg.sigma_perm.(fld).idx(:, idx_pcor(i)))-1), [Ax.YLim(1), Ax.YLim(1)].*0.5, '-k', 'LineWidth', 3);
    end
end

if strcmpi(fld, 'aw')
    YData = double(perms(1).y.pbo.pval < 0.05);
    YData(~YData) = nan;
    plot(Ax, XData, YData*(Ax.YLim(1).*0.33), '-', 'LineWidth', 3, 'Color', css_standard_colors('pbo'))

    YData = double(perms(1).y.etc.pval < 0.05);
    YData(~YData) = nan;
    plot(Ax, XData, YData*(Ax.YLim(1).*0.67), '-', 'LineWidth', 3, 'Color', css_standard_colors('etc'))
end

% For placebo
[~, YData] = withinSubMean(T(pcfg.idx.pbo.(fld), :), 'smtdata', 1);
YData = YData';
EData = tinv(0.975, size(YData, 2)-1).*(std(YData, [], 2)./sqrt(size(YData, 2)));
errorpatch(Ax, XData, mean(YData, 2), EData, ...
    'FaceColor', css_standard_colors('pbo'), ...
    'FaceAlpha', 0.3);
p(1) = plot(Ax, XData, mean(YData, 2), '-', ...
    'LineWidth', 1, ...
    'Color', css_standard_colors('pbo'));

% For THC/CBD
[~, YData] = withinSubMean(T(pcfg.idx.etc.(fld), :), 'smtdata', 1);
YData = YData';
EData = tinv(0.975, size(YData, 2)-1).*(std(YData, [], 2)./sqrt(size(YData, 2)));
errorpatch(Ax, XData, mean(YData, 2), EData, ...
    'FaceColor', css_standard_colors('etc'), ...
    'FaceAlpha', 0.3);
p(2) = plot(Ax, XData, mean(YData, 2), '-', ...
    'LineWidth', 1, ...
    'Color', css_standard_colors('etc'));

if strcmpi(fld, 'cs')
    l = legend(p, {'placebo', 'THC/CBD'}, ...
        'Box', 'off', ...
        'FontSize', 8, ...
        'Location', 'eastoutside');
    l.Position(1:2) = [-0.0175, 0.4];
    l.ItemTokenSize = [10,10];
end

end