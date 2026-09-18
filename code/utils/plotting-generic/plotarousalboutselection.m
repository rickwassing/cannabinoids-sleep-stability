function Fig = plotarousalboutselection(EEG)

Fig = figure('Position', [1, 860, 1920, 130]);
Ax = axes(Fig);
Ax.NextPlot = 'add';
Ax.TickLength = [0.001, 0.001];
for i = 1:length(EEG.event)
    if ~EEG.event(i).is_selected
        continue
    end
    if EEG.event(i).is_awakening
        clr = standard_colors('red');
    else
        clr = [0.3, 0.3, 0.3];
    end
    plot(Ax, [EEG.event(i).latency, EEG.event(i).latency+EEG.event(i).duration]./(24*60*60*EEG.srate), [1.5, 1.5], ...
        'Marker', 'none', ...
        'Color', clr, ...
        'LineWidth', 2)
    plot(Ax, [EEG.event(i).latency, EEG.event(i).latency]./(24*60*60*EEG.srate), [-4.5, 1.5], ...
        'Marker', 'none', ...
        'LineStyle', '-', ...
        'Color', clr.^0.25, ...
        'LineWidth', 0.5)
end
plotHypnogram(Ax, EEG);
end