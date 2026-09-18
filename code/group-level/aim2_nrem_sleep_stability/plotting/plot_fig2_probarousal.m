function [Ax] = plot_fig2_probarousal(Fig, T, perms, fld, pcfg)

% Create axes
Ax = axes(Fig, 'NextPlot', 'add');
Ax.Box = 'on';
Ax.FontSize = 8;
Ax.TickLength = [0 0];
Ax.YLim = [0 0.1];
Ax.YTick = [0 0.1];
Ax.YTickLabel{1} = '';
Ax.XTickLabel = {};
Ax.XGrid = 'on';
Ax.Color = [0.96 0.97 0.99];

switch fld
    case 'aw'
        Ax.Title.String = 'AWAKENINGS';
        Ax.Title.Color = css_standard_colors('aw');
        Ax.YLabel.String = 'Pr(aro)';
        Ax.YLabel.FontSize = 8;
    case 'cs'
        Ax.Title.String = 'CONT. SLEEP';
        Ax.Title.Color = css_standard_colors('cs');
        Ax.YLabel.String = ' ';
        Ax.YLabel.FontSize = 8;
        Ax.YTickLabel{2} = '';
end
Ax.Title.FontWeight = 'normal';
Ax.Title.FontSize = 10;

% Time series
XData = (-59.9:0.1:0)-(abs(pcfg.const.crop_aro(2))/10-30);

% Plot arousal marker
plot(Ax, [0, 0], Ax.YLim, ':k')

% Plot significant bins
YData = double(perms(1).pr.(fld).pval < 0.05);
YData(~YData) = nan;
plot(Ax, XData, YData*(Ax.YLim(2).*0.9), '-', 'LineWidth', 3, 'Color', css_standard_colors('bluegrey'))

% Plot significant bins
if strcmpi(fld, 'aw')
    YData = double(perms(1).pr.pbo.pval < 0.05);
    YData(~YData) = nan;
    plot(Ax, XData, YData*(Ax.YLim(2).*0.9), '-k', 'LineWidth', 3, 'Color', css_standard_colors('pbo'))
end
if strcmpi(fld, 'cs')
    YData = double(perms(1).pr.etc.pval < 0.05);
    YData(~YData) = nan;
    plot(Ax, XData, YData*(Ax.YLim(2).*0.9), '-k', 'LineWidth', 3, 'Color', css_standard_colors('etc'))
end

% Plot Pr(Aro) for placebo
[PrAroData] = mean(withinSubMean(T(pcfg.idx.pbo.(fld), :), 'pr_aro'));
patch(Ax, ...
    'XData', [XData(1), XData, XData(end)], ...
    'YData', [0, PrAroData, 0], ...
    'LineStyle', 'none', ...
    'FaceColor', css_standard_colors('pbo'), ...
    'FaceAlpha', 0.3)

% Plot Pr(Aro) for THC/CBD
[PrAroData] = mean(withinSubMean(T(pcfg.idx.etc.(fld), :), 'pr_aro'));
patch(Ax, ...
    'XData', [XData(1), XData, XData(end)], ...
    'YData', [0, PrAroData, 0], ...
    'LineStyle', 'none', ...
    'FaceColor', css_standard_colors('etc'), ...
    'FaceAlpha', 0.3)

end