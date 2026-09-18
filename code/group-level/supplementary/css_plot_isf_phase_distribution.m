function css_plot_isf_phase_distribution()

load('analysis_2a.mat', 'ANGA');

XData = 0:pi/180:100*pi;
YData = sin(XData);
AData = angle(hilbert(YData));

Fig = figure('Units','centimeters');
Fig.Position = [15 15 8 6];


Ax = axes(Fig, 'NextPlot', 'add');

Ax.Box = 'on';
Ax.FontSize = 8;
Ax.Layer = 'top';
Ax.XLim = [48*pi, 52*pi] + [-0.08*pi 0.08*pi];
Ax.YLim = [-pi, pi] + [-0.08*pi 0.08*pi];
Ax.XTick = (-3*pi:0.5*pi:3*pi);
Ax.YTick = [-pi, -0.5*pi, 0, 0.5*pi, pi];
Ax.YTickLabel = {'-\pi', '-0.5\pi','0', '0.5\pi', '\pi'};
Ax.Color = [0.95, 0.96, 0.98];
Ax.TickLength = [0, 0];
Ax.XLabel.String = 'time';
Ax.XLabel.FontSize = 10;
Ax.YLabel.String = ' sin(t)';
Ax.YLabel.FontSize = 10;
Ax.OuterPosition = [0 0.5 1 0.5];

plot(Ax, [0, 314], [0 0], '-k', 'LineWidth', 0.25);
plot(Ax, XData, YData, '-', 'Color', [0.839, 0.847, 0.863], ...
    'LineWidth', 3)
plot(Ax, XData, AData, ':', 'Color', 'k', ...
    'LineWidth', 1)

pcfg = struct();
pcfg.resolution = 30;
pcfg.nbins = 360/pcfg.resolution;
pcfg.bins.edges = linspace(-pi, pi, pcfg.nbins+1);
pcfg.bins.centers = pcfg.bins.edges(1:end-1) + diff(pcfg.bins.edges)/2;

ANGA = sort(ANGA(:));

close all

Fig = figure('Units','centimeters');
Fig.Position = [15 15 8 8];

clear Ax

% Create axes object
Ax(1) = axes(Fig, 'NextPlot', 'add');

% Set axis properties
Ax(1).Box = 'on';
Ax(1).FontSize = 8;
Ax(1).Layer = 'top';
Ax(1).XLim = [0, 2*pi] + [-0.19*pi 0.19*pi];
Ax(1).YLim = [-1.08*pi, 1.08*pi];
Ax(1).XTick = (-3*pi:0.5*pi:3*pi);
Ax(1).XGrid = 'on';
Ax(1).YTick = [-pi, 0, pi];
Ax(1).Color = [0.95, 0.96, 0.98];
Ax(1).TickLength = [0, 0];
Ax(1).YLabel.String = ' ';
Ax(1).YLabel.FontSize = 10;
Ax(1).XTickLabel = {'180', '270', '0', '90', '180', '270', '0', '90', '180', '270', '0', '90', '180'};
Ax(1).YTickLabel = {'-\pi', '0', '\pi'};
Ax(1).OuterPosition = [0, 0.7, 1, 0.3];

plot(Ax(1), [-3*pi 3*pi], [0, 0], '-k', 'LineWidth', 0.25)
plot(Ax(1), -3*pi:pi/180:3*pi, cos(-3*pi:pi/180:3*pi).*1.5, '-k', 'LineWidth', 3)
plot(Ax(1), -3*pi:pi/180:3*pi, angle(hilbert(cos(-3*pi:pi/180:3*pi))), ':k', 'LineWidth', 1.5)

% Create axes object
Ax(2) = axes(Fig, 'NextPlot', 'add');

% Set axis properties
Ax(2).Box = 'on';
Ax(2).FontSize = 8;
Ax(2).Layer = 'top';
Ax(2).XGrid = 'on';
Ax(2).XLim = [0, 2*pi] + [-0.19*pi 0.19*pi];
Ax(2).XTick = (-3*pi:0.5*pi:3*pi);
Ax(2).Color = [0.95, 0.96, 0.98];
Ax(2).TickLength = [0, 0];
Ax(2).XTickLabel = {'trough', 'asc', 'peak', 'desc', 'trough', 'asc', 'peak', 'desc', 'trough', 'asc', 'peak', 'desc', 'trough'};
Ax(2).XLabel.String = '\sigma ISF Phase';
Ax(2).XLabel.FontSize = 10;
Ax(2).YLabel.String = ' PDF';
Ax(2).YLabel.FontSize = 10;
Ax(2).OuterPosition = [0, 0, 1, 0.7];

% Extract the phase angle 
AData = ANGA(:);
AData(AData > max(pcfg.bins.edges)) = AData(AData > max(pcfg.bins.edges)) - 2*pi;
Pr = histcounts(AData, 'BinEdges', pcfg.bins.edges, 'Normalization', 'pdf');

% Plot bar graphs and set their face-color
XData = pcfg.bins.centers;
XData = [XData-2*pi, XData, XData+2*pi]; 
YData = Pr';
YData = [YData; YData; YData];

h = bar(XData, YData, ...
    'BarWidth', 0.8, ...
    'GroupWidth', 1, ...
    'LineStyle', 'none');

h(1).FaceColor = standard_colors('bluegrey');

Ax(1).Position([1 3]) = Ax(2).Position([1 3]);

exportgraphics(Fig, './figures/supp_phaseangledistribution.png', 'Resolution', 600)

end