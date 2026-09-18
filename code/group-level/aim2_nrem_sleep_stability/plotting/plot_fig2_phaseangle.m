function [Ax] = plot_fig2_phaseangle(Fig, T, fld, delay, pcfg)

doAvg = true;

% Create axes object
Ax = axes(Fig, 'NextPlot', 'add');

% Set axis properties
Ax.Box = 'on';
Ax.FontSize = 8;
Ax.Layer = 'top';
Ax.XLim = [0, 2*pi] + [-0.08*pi 0.08*pi];
Ax.YLim = pcfg.ylim_phase;
Ax.YTick = [0, Ax.YLim(2)];
Ax.XTick = (-3*pi:0.5*pi:3*pi);
Ax.Color = [0.95, 0.96, 0.98];
Ax.TickLength = [0, 0];
Ax.XTickLabel = {'trough', 'asc', 'peak', 'desc', 'trough', 'asc', 'peak', 'desc', 'trough', 'asc', 'peak', 'desc', 'trough'};
Ax.XLabel.String = '\sigma ISF Phase';
Ax.XLabel.FontSize = 10;
switch fld
    case 'aw'
        Ax.YLabel.String = ' # arousals';
        Ax.YLabel.FontSize = 10;
        Ax.Title.String = 'AWAKENINGS';
        Ax.Title.Color = standard_colors('aw');
    case 'cs'
        Ax.YLabel.String = ' ';
        Ax.YLabel.FontSize = 10;
        Ax.YTick = [];
        Ax.Title.String = 'CONT. SLEEP';
        Ax.Title.Color = standard_colors('cs');
end
Ax.Title.FontWeight = 'normal';
Ax.Title.FontSize = 10;

% plot XGrid
XData = [...
    -3.0*pi, -3.0*pi, nan, ...
    -2.5*pi, -2.5*pi, nan, ...
    -2.0*pi, -2.0*pi, nan, ...
    -1.5*pi, -1.5*pi, nan, ...
    -1.0*pi, -1.0*pi, nan, ...
    -0.5*pi, -0.5*pi, nan, ...
    0.0*pi, 0.0*pi, nan, ...
    0.5*pi, 0.5*pi, nan, ...
    1.0*pi, 1.0*pi, nan, ...
    1.5*pi, 1.5*pi, nan, ...
    2.0*pi, 2.0*pi, nan, ...
    2.5*pi, 2.5*pi, nan, ...
    3.0*pi, 3.0*pi, nan, ...
    ];
YData = repmat([0, Ax.YLim(2), nan], 1, 13);
plot(Ax, XData+pi/4, YData, '-', 'Color', [0.839, 0.847, 0.863])

% Extract the phase angle of arousals leading to awakenings

% - Placebo
if doAvg
    AData_pbo = within_chan_circ_mean(T(pcfg.idx.pbo.(fld), :), delay);
else
    AData_pbo = T.(delay)(pcfg.idx.pbo.(fld), :);
end
AData_pbo(AData_pbo > max(pcfg.bins.edges)) = AData_pbo(AData_pbo > max(pcfg.bins.edges)) - 2*pi;
Pr_pbo = histcounts(AData_pbo, 'BinEdges', pcfg.bins.edges);

% - THC/CBD
if doAvg
    AData_etc = within_chan_circ_mean(T(pcfg.idx.etc.(fld), :), delay);
else
    AData_etc = T.(delay)(pcfg.idx.etc.(fld), :);
end
AData_etc(AData_etc > max(pcfg.bins.edges)) = AData_etc(AData_etc > max(pcfg.bins.edges)) - 2*pi;
Pr_etc = histcounts(AData_etc, 'BinEdges', pcfg.bins.edges);

% Plot ISF phase angle
plot(Ax, -3*pi:pi/72:3*pi, 0.5*Ax.YLim(2)+cos((-3*pi:pi/72:3*pi))*0.25*Ax.YLim(2),  '-', ...
    'Color', [0.85, 0.86, 0.88], ...
    'LineWidth', 3)

% Plot bar graphs and set their face-color
XData = pcfg.bins.centers;
XData = [XData-2*pi, XData, XData+2*pi]; 
YData = [Pr_pbo', Pr_etc'];
YData = [YData; YData; YData];

h = bar(XData, YData, ...
    'BarWidth', 0.8, ...
    'GroupWidth', 1, ...
    'LineStyle', 'none');

h(1).FaceColor = standard_colors('pbo');
h(2).FaceColor = standard_colors('etc');

end