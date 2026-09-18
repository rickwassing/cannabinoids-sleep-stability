function plotisfspect(ISF)

YData = [ISF.data, nan(size(ISF.data, 1), 1)]';
YData = YData(:)';
XData = repmat([asrow(ISF.freqs), nan], 1, size(ISF.data, 1));

close all 
figure
ax = axes('NextPlot','add');
patch('XData', XData, 'YData', YData, 'EdgeColor', 'k', 'EdgeAlpha', 0.05)
plot(ISF.features(1).data, ISF.features(2).data, '.r', 'MarkerSize', 5)
ax.XLim = [0, 0.1];

end