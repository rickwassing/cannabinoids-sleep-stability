function [Ax] = plot_fig2_sigmatrace(Fig, SIG, FSIG, fld, pcfg)

% Get indexex
idx_ar = round(SIG.event(pcfg.idx.aro.(fld)).latency);
idx_x = idx_ar+pcfg.xlim(1)*SIG.srate:idx_ar+pcfg.xlim(2)*SIG.srate;

% Create axes
Ax = axes(Fig, 'NextPlot', 'add');
Ax.XAxis.Visible = 'off';
Ax.YAxis.Color = 'w';
Ax.XLim = [idx_x(1), idx_x(end)]./SIG.srate;
Ax.YLim = pcfg.ylim_sigma;
Ax.YTick = [];
Ax.YLabel.String = '\sigma';
Ax.YLabel.FontSize = 12;
Ax.YLabel.Color = 'k';

% Plot arousal marker
XData = SIG.times(idx_ar);
plot(Ax, [XData, XData], pcfg.ylim_sigma, ':k')

% Extract data for plotting
null_r = pcfg.xlim(2)*SIG.srate+(abs(pcfg.const.crop_aro(2))-30*SIG.srate);
null_f = pcfg.xlim(2)*SIG.srate;
XData = SIG.times(idx_x);
XData_r = XData;
XData_f = XData;
XData_f(end-null_f+1:end) = nan;
YData_r = zscore(SIG.data(idx_x));
YData_f = zscore(FSIG.data(idx_x))*2+2;
YData_a = angle(hilbert(FSIG.data));
YData_a = YData_a(idx_x);

% Plot sigma power, filtered data and angle
clear p l
p(1) = plot(Ax, XData_r, YData_r, '-', 'Color', [0.65 0.66 0.68], 'LineWidth', 1);
patch(Ax, 'XData', [XData, XData(end)+1, -1], 'YData', [YData_r, -3.3, -3.3], 'CData', [YData_r, -3.3, -3.3], ...
    'EdgeColor', standard_colors('bluegrey'), ...
    'FaceColor', 'interp');
p(3) = plot(Ax, XData_f, YData_a, '-', 'Color', [0 0 0], 'LineWidth', 1.5);
p(2) = plot(Ax, XData_f, YData_f, ':', 'Color', [0 0 0], 'LineWidth', 1.5);
plot(Ax, max(XData_f), YData_a(find(~isnan(XData_f), 1, 'last')), '.', 'Color', [0 0 0], 'MarkerSize', 12);

if strcmpi(fld, 'aw')
    % Plot time scale
    plot(Ax, [XData(1)+55, XData(1)+65], [5, 5], '-k', 'LineWidth', 2)
    text(Ax, XData(1)+60, 5, '10 s', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 8)
end

if strcmpi(fld, 'cs')
    % Legend
    l = legend(p, {'raw', 'filtered', 'phase'}, 'FontSize', 8, 'Box', 'off');
    l.ItemTokenSize = [10 10];
    l.Position(1:2) = [0.88, 0.71];
end

CMap = standard_colors('bluegrey');
Ax.UserData.CMap = [linspace(1, CMap(1), 256); linspace(1, CMap(2), 256); linspace(1, CMap(3), 256)]';

end