function plotallsigmachannels(SIGMA)

YScale = median(SIGMA.data(:), 'omitnan') + 24*mad(SIGMA.data(:));
boundaries = [SIGMA.event(strcmpi({SIGMA.event.type}, 'boundary')).latency] ./ SIGMA.srate;

close all
Fig = figure(...
    'Position', [1, 62, 1920, 915]);

Ax = axes(Fig, ...
    'XLim', [0, SIGMA.xmax], ...
    'YLim', [-0.5, 45.5], ...
    'TickLength', [0, 0], ...
    'Layer', 'top', ...
    'Box', 'on', ...
    'NextPlot', 'add', ...
    'TickDir', 'out', ...
    'OuterPosition', [0, 0, 1/4, 1]+[-0.025 0 0.05 0]);
Ax.Title.String = sprintf('Y-Scale: %.1f uV^2', YScale);

for b = 1:length(boundaries)
    plot([boundaries(b), boundaries(b)], [-0.5 45.5], '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.5)
end

j = 1;
cnt = 0;
col = 0;
for i = 1:SIGMA.nbchan
    if mod(i, 45) == 0
        cnt = 0;
        col = col+1;
        Ax.YTick = 0:length(j:i-1)-1;
        Ax.YTickLabel = {SIGMA.chanlocs(j:i-1).labels};
        Ax = axes(Fig, ...
            'XLim', [0, SIGMA.xmax], ...
            'YLim', [-0.5, 45.5], ...
            'TickLength', [0, 0], ...
            'Layer', 'top', ...
            'Box', 'on', ...
            'NextPlot', 'add', ...
            'TickDir', 'out', ...
            'OuterPosition', [col/4, 0, 1/4, 1]+[-0.025 0 0.05 0]);
        for b = 1:length(boundaries)
            plot([boundaries(b), boundaries(b)], [-0.5 45.5], '-', 'Color', [0.75 0.75 0.75], 'LineWidth', 0.5)
        end
        j = i;
    end
    cnt = cnt+1;
    XData = SIGMA.times;
    YData = SIGMA.data(i, :);
    EData = nan(1, length(YData));
    EData(YData > YScale) = cnt-1;
    MinYData = min(YData);
    [MaxYData, MaxXDataIdx] = max(YData);
    YData = YData ./ YScale;
    YData = YData + cnt-1;
    plot(XData, YData, '-k');
    plot(XData, EData, '.r', 'MarkerSize', 10)
end
Ax.YTick = 0:length(j:i)-1;
Ax.YTickLabel = {SIGMA.chanlocs(j:i).labels};

exportgraphics(Fig, sprintf('figures/subject-level/sigmapowertimeseries/%s.png', SIGMA.setname), 'Resolution', 600)


