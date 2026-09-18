Files = dir('derivatives/EEG-preproc/sub-*/ses-*/sub-*-*-psg_desc-preprochr_hr.set');

for i = 33
    HR = LoadDataset(fullfile(Files(i).folder, Files(i).name), 'all');
    hyp = css_eeglab2hypnogram(HR);

    close all
    clear Ax
    Fig = figure('Position', [1 680 1920 260]);
    Ax(1) = axes(Fig, ...
        'OuterPosition', [0 0.5 1 0.5], ...
        'TickLength', [0, 0]);
    plotHypnogram(Ax(1), HR)

    bouts = getnrembouts(hyp, HR.srate, 300);
    for j = 1:size(bouts)
        YData = [0.875; 0.875; 1.125; 1.125] - 3.5;
        XData = [bouts(j, 1); bouts(j, 2); bouts(j, 2); bouts(j, 1)] ./ (24*60*60);
        Vertices = [XData, YData];
        Faces = 1:4;
        CData = [standard_colors('orange'); standard_colors('orange'); standard_colors('orange'); standard_colors('orange')];
        patch(Ax(1), 'Faces', Faces, 'Vertices', Vertices, ...
            'FaceVertexCData', CData, ...
            'FaceColor','interp', ...
            'LineStyle', 'none');
    end


    Ax(2) = axes(Fig, ...
        'OuterPosition', [0 0 1 0.5], ...
        'TickLength', [0, 0], ...
        'NextPlot', 'add');

    plot(Ax(2), HR.times./(24*60*60), HR.data, '-k')

    for j = 1:size(bouts)
        YData = [40; 40; 120; 120];
        XData = [bouts(j, 1); bouts(j, 2); bouts(j, 2); bouts(j, 1)] ./ (24*60*60);
        Vertices = [XData, YData];
        Faces = 1:4;
        CData = [standard_colors('orange'); standard_colors('orange'); standard_colors('orange'); standard_colors('orange')];
        patch(Ax(2), 'Faces', Faces, 'Vertices', Vertices, ...
            'FaceVertexCData', CData, ...
            'FaceAlpha', 0.4, ...
            'FaceColor','interp', ...
            'LineStyle', 'none');
    end


    Ax(1).XLim = [HR.times(1), HR.times(end)]./(24*60*60);
    Ax(2).XLim = [HR.times(1), HR.times(end)]./(24*60*60);
    Ax(2).XTick = Ax(1).XTick;
    Ax(2).XTickLabel = Ax(1).XTickLabel;
    Ax(2).TickLength = [0 0];
    Ax(2).YLim = [0 120];
    Ax(2).YTick = 40:10:120;
    Ax(1).XGrid = 'on';
    Ax(1).Color = [0.96 0.97 0.99];
    Ax(2).XGrid = 'on';
    Ax(2).YGrid = 'on';
    Ax(2).Color = [0.96 0.97 0.99];

    linkaxes(Ax, 'x')

    exportgraphics(Fig, ['inspect/ecg/', HR.setname, '.png'], 'Resolution', 300)
end