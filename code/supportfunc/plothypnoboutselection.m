function Fig = plothypnoboutselection(EEG, bouts)
close all;
Fig = figure('Position', [1, 860, 1920, 130]);
Ax = axes(Fig);
Ax.TickLength = [0.001, 0.001];
for i = 1:size(bouts)
    % YData = [0.875; 0.875; 1.125; 1.125] - 3.5;
    YData = [-5.5; -5.5; 1.5; 1.5];
    XData = [bouts(i, 1); bouts(i, 2); bouts(i, 2); bouts(i, 1)] ./ (24*60*60);
    Vertices = [XData, YData];
    Faces = 1:4;
    CData = [standard_colors('orange'); standard_colors('orange'); standard_colors('orange'); standard_colors('orange')];
    patch(Ax, 'Faces', Faces, 'Vertices', Vertices, ...
        'FaceVertexCData', CData, ...
        'FaceColor', 'interp', ...
        'FaceAlpha', 0.1, ...
        'LineStyle', 'none');
end
plotHypnogram(Ax, EEG);
end