function Ax = plot_fig3_sigmatrace(Fig, H, pcfg)

idx.sel = 109;
idx.chan = 2;
idx.nrem.ev = H(idx.sel).event.channel == idx.chan & strcmpi(H(idx.sel).event.stage, 'nrem');
idx.rem.ev = H(idx.sel).event.channel == idx.chan & strcmpi(H(idx.sel).event.stage, 'rem');
idx.nrem.trough = H(idx.sel).event.latency(idx.nrem.ev);
idx.nrem.peak = H(idx.sel).event.peak_latency(idx.nrem.ev);
idx.rem.trough = H(idx.sel).event.latency(idx.rem.ev);
idx.rem.peak = H(idx.sel).event.peak_latency(idx.rem.ev);

Ax = axes(Fig, 'NextPlot', 'add');
Ax.XAxis.Visible = 'off';
Ax.YAxis.Color = 'w';
Ax.XLim = [H(idx.sel).times(1), H(idx.sel).times(end)];
Ax.YTick = [];
Ax.YLim = [-2.5 2.5];
Ax.YLabel.String = '\sigma';
Ax.YLabel.FontSize = 12;
Ax.YLabel.Color = 'k';

% Plot sigma
XData = asrow(H(idx.sel).times);
YData = asrow(H(idx.sel).rawsigma(:, idx.chan));
YData = 10.^(YData ./ 10)-2;
p(1) = plot(Ax, H(idx.sel).times, YData, '-', 'Color', [1, 1, 1], 'LineWidth', 0.25);
patch(Ax, 'XData', [XData, XData(end)+1, XData(1)-1], 'YData', [YData, -2.53, -2.53], 'CData', [YData, -3.3, -3.3], ...
    'EdgeColor', css_standard_colors('bluegrey').^3.0, ...
    'FaceColor', 'interp');
p(2) = plot(Ax, H(idx.sel).times, H(idx.sel).filtsigma(:, idx.chan), ':', 'Color', [0, 0, 0], 'LineWidth', 1.5);

% Plot REM onset
plot(Ax, [6, 6], Ax.YLim, ':k', 'LineWidth', 1)

for j = 1:length(idx.nrem.trough)
    plot(Ax, ...
        [H(idx.sel).times(idx.nrem.trough(j)), H(idx.sel).times(idx.nrem.trough(j)), H(idx.sel).times(idx.nrem.peak(j))], ...
        [H(idx.sel).filtsigma(idx.nrem.trough(j), idx.chan), H(idx.sel).filtsigma(idx.nrem.peak(j), idx.chan), H(idx.sel).filtsigma(idx.nrem.peak(j), idx.chan)], ...
        'LineStyle', '-', ...
        'Color', standard_colors('blue'))
end

plot(Ax, ...
    [H(idx.sel).times(idx.rem.trough(1)), H(idx.sel).times(idx.rem.trough(1)), H(idx.sel).times(idx.rem.peak(1))], ...
    [H(idx.sel).filtsigma(idx.rem.trough(1), idx.chan), H(idx.sel).filtsigma(idx.rem.peak(1), idx.chan), H(idx.sel).filtsigma(idx.rem.peak(1), idx.chan)], ...
    'LineStyle', '-', ...
    'Color', standard_colors('red'))

% Plot trough and peak markers
plot(Ax, H(idx.sel).times(idx.nrem.trough), H(idx.sel).filtsigma(idx.nrem.trough, idx.chan), 'o', 'MarkerSize', pcfg.markersize, 'Color', standard_colors('blue'), 'MarkerFaceColor', standard_colors('blue'), 'MarkerEdgeColor', 'w')
plot(Ax, H(idx.sel).times(idx.nrem.peak), H(idx.sel).filtsigma(idx.nrem.peak, idx.chan), 'o', 'MarkerSize', pcfg.markersize, 'Color', standard_colors('blue'), 'MarkerFaceColor', standard_colors('blue'), 'MarkerEdgeColor', 'w')
plot(Ax, H(idx.sel).times(idx.rem.trough), H(idx.sel).filtsigma(idx.rem.trough, idx.chan), 'o', 'MarkerSize', pcfg.markersize, 'Color', standard_colors('red'), 'MarkerFaceColor', standard_colors('red'), 'MarkerEdgeColor', 'w')
plot(Ax, H(idx.sel).times(idx.rem.peak), H(idx.sel).filtsigma(idx.rem.peak, idx.chan), 'o', 'MarkerSize', pcfg.markersize, 'Color', standard_colors('red'), 'MarkerFaceColor', standard_colors('red'), 'MarkerEdgeColor', 'w')

% Plot time scale
plot(Ax, [10 60], [Ax.YLim(1) Ax.YLim(1)], '-k', 'LineWidth', 2)
text(Ax, 35, Ax.YLim(1), '50 s', 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 8)

% Plot text
text(Ax, 6, Ax.YLim(2), ' N2 → REM ', ...
    'FontSize', 8, ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'top', ...
    'Color', 'k')

CMap = standard_colors('bluegrey');
Ax.Colormap = [linspace(1, CMap(1), 256); linspace(1, CMap(2), 256); linspace(1, CMap(3), 256)]';
Ax.UserData.CMap = [linspace(1, CMap(1), 256); linspace(1, CMap(2), 256); linspace(1, CMap(3), 256)]';

end