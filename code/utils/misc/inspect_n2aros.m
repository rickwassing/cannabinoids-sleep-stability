% NOTE: no current call sites found in code/ as of 2026-09-18 (Phase 3
% refactor). This is a bare script (not a function), likely an
% interactive/manual QC snippet. Kept in utils/misc/ rather than archived,
% since unclear current usage does not necessarily mean unused; flagged
% for your review.

ARO = [];

Files = dir('./derivatives/EEG-preproc/sub-*/ses-*/sub-*_desc-preproc_eeg.set');

for f = 30:length(Files)
    EEG = LoadDataset(fullfile(Files(f).folder, Files(f).name), 'all');

    cfg = struct();
    cfg.stage = -2;
    cfg.cutoff = 100;
    cfg.allowalpha = false;
    cfg.allowmicro = false;
    % Get hypnogram table
    HYP = css_eeglab2hypnogram(EEG);
    % Find the arousal bouts
    [EEG, arousalbouts, ids] = getarousalbouts(EEG, HYP, cfg);
    arousalbouts = round(arousalbouts);

    close all

    Fig = figure('Color', 'w', 'Position', [0 150 1920 250]);

    clear Ax
    Ax(1) = axes('NextPlot', 'add');
    Ax(1).OuterPosition = [0 0.67 1 0.33];
    Ax(1).TickLength = [0 0];
    hyp = plotHypnogram(Ax(1), EEG);
    Ax(1).XTick = [];
    Ax(1).XTickLabel = {''};
    Ax(1).XGrid = 'on';
    Ax(1).YLim = [-3.5 1.75];

    Ax(2) = axes('NextPlot', 'add');
    Ax(2).OuterPosition = [0 0 1 0.67];
    Ax(2).TickLength = [0 0];
    Ax(2).XTick = [];
    Ax(2).YLim = [-150 150];
    XData = [arousalbouts(:, 2)-30*EEG.srate, arousalbouts(:, 2)-30*EEG.srate, nan(size(arousalbouts, 1), 1)]'./EEG.srate;
    YData = [ones(1, size(XData, 2)); -ones(1, size(XData, 2)); nan(1, size(XData, 2))];
    plot(Ax(2), XData(:), YData(:).*150, ':k', 'LineWidth', 1)
    plot(Ax(2), EEG.times, EEG.data(94, :), '-k', 'LineWidth', 0.25)

    for i = 1:size(arousalbouts, 1)

        try

            Ax(1).Title.String = sprintf('%s ID %i', strrep(EEG.setname, '_', ' '), ids(i));
            Ax(1).Title.FontSize = 10;
            Ax(1).Title.FontWeight = 'normal';
            Ax(1).XLim = [(EEG.times(arousalbouts(i, 2))-90)/(60*60*24), (EEG.times(arousalbouts(i, 2))+30)/(60*60*24)];
            Ax(1).XTick = [(EEG.times(arousalbouts(i, 2))-45)/(60*60*24), (EEG.times(arousalbouts(i, 2))-15)/(60*60*24)];
            Ax(2).XLim = [EEG.times(arousalbouts(i, 2))-45, EEG.times(arousalbouts(i, 2))-15];

            aro = table();
            aro.setname = {EEG.setname};
            aro.id = ids(i);
            aro.include = 1;

            exportgraphics(Fig, sprintf('./inspect/n2aros/%s_id-%i.png', EEG.setname, ids(i)), 'Resolution', 144)

            ARO = [ARO; aro];
        catch
            continue
        end

    end
end

writetable(ARO, './inspect/n2aros.csv')