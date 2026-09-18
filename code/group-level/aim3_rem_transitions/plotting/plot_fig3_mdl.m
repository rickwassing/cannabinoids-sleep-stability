function Ax = plot_fig3_mdl(Fig, desmat, fld)

Ax = axes(Fig, 'NextPlot', 'add');
if strcmpi(fld, 'duration')
    h.a = plot(Ax, fitlm(desmat(strcmpi(desmat.stage, 'rem'), :), sprintf('remeplat ~ 1 + %s', fld)));
    h.b = plot(Ax, fitlm(desmat(strcmpi(desmat.stage, 'nrem'), :), sprintf('remeplat ~ 1 + %s', fld)));
else
    h.a = plot(Ax, fitlm(desmat, sprintf('remeplat ~ 1 + %s', fld)));
end
x = scatter(Ax, h.a(1).XData, h.a(1).YData, ...
    'SizeData', 5, ...
    'Marker', 'o', ...
    'MarkerEdgeColor', 'none', ...
    'MarkerFaceColor', 'k', ...
    'MarkerFaceAlpha', 0.1);
if isfield(h, 'b')
    x = scatter(Ax, h.b(1).XData, h.b(1).YData, ...
    'SizeData', 5, ...
    'Marker', 'o', ...
    'MarkerEdgeColor', 'none', ...
    'MarkerFaceColor', 'k', ...
    'MarkerFaceAlpha', 0.1);
end


if isfield(h, 'b')
    h.a(2).Color = standard_colors('red');
    h.a(3).Color = standard_colors('red');
    h.b(2).Color = standard_colors('blue');
    h.b(3).Color = standard_colors('blue');
else
    h.a(2).Color = 'k';
    h.a(3).Color = 'k';
end
try
    h.a(4).Color = 'k';
    if isfield(h, 'b')
        h.b(4).Color = 'k';
    end
end
h.a(2).LineWidth = 1;
h.a(3).LineWidth = 0.25;
if isfield(h, 'b')
    h.b(2).LineWidth = 1;
    h.b(3).LineWidth = 0.25;
end
try
    h.a(4).LineWidth = 0.25;
    if isfield(h, 'b')
        h.b(4).LineWidth = 0.25;
    end
end
h.a(3).LineStyle = '-';
if isfield(h, 'b')
    h.b(3).LineStyle = '-';
end
try
    h.a(4).LineStyle = '-';
    if isfield(h, 'b')
        h.b(4).LineStyle = '-';
    end
end
delete(Ax.Legend)
switch fld
    case 'amplitude'
        Ax.XLabel.String = 'Amplitude (dB)';
        Ax.XLim = [0, 3];
    case 'duration'
        Ax.XLabel.String = 'Half-period (s)';
        Ax.XLim = [5, 65];
        Ax.XTick = 0:10:60;
        Ax.XTickLabelRotation = 0;
end
Ax.XLabel.FontSize = 10;
Ax.YLabel.String = 'REML (min)';
Ax.YLabel.FontSize = 10;
Ax.Title.String = '';
Ax.Color = [0.96 0.97 0.99];
Ax.Box = 'on';
Ax.FontSize = 8;
Ax.TickLength = [0, 0];
Ax.YTick = 0:60:360;
Ax.YLim = [0, 420];
Ax.YGrid = 'on';

children = get(Ax, 'Children');
set(Ax, 'Children', flipud(children));


delete(h.a(1))
if isfield(h, 'b')
    delete(h.b(1))
end

end