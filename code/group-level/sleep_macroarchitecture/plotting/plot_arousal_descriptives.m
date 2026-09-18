function plot_arousal_descriptives(PSG, ARO, S)
% -------------------------------------------------------------------------
% Figure: Arousal descriptives (arousal index by stage, duration histogram,
% and probability-of-awakening-by-duration for N1/N2/N3).
% -------------------------------------------------------------------------
BinEdges   = 1.5:3:40.5;
BinCenters = BinEdges(1:end-1) + mean(diff(BinEdges))/2;

close all
Fig = figure('Units','centimeters','Position',[5 5 8.5 8], 'Color', 'w');

clear Ax; ai = 0;

% Panel A — Arousal Index by Stage
ai = ai + 1;
Ax(ai) = axes(Fig, 'NextPlot', 'add');

idxPBO = strcmpi(PSG.condition, 'placebo');

YData = [ ...
    mean(PSG.i_aro_n1(idxPBO)), mean(PSG.i_aro_n1(~idxPBO)); ...
    mean(PSG.i_aro_n2(idxPBO)), mean(PSG.i_aro_n2(~idxPBO)); ...
    mean(PSG.i_aro_n3(idxPBO)), mean(PSG.i_aro_n3(~idxPBO))];

EData = tinv(0.975, length(idxPBO)-1) .* [ ...
    std(PSG.i_aro_n1(idxPBO))/sqrt(sum(idxPBO)), std(PSG.i_aro_n1(~idxPBO))/sqrt(sum(~idxPBO)); ...
    std(PSG.i_aro_n2(idxPBO))/sqrt(sum(idxPBO)), std(PSG.i_aro_n2(~idxPBO))/sqrt(sum(~idxPBO)); ...
    std(PSG.i_aro_n3(idxPBO))/sqrt(sum(idxPBO)), std(PSG.i_aro_n3(~idxPBO))/sqrt(sum(~idxPBO))];

b = bar(Ax(ai), 1:size(YData,1), YData);
b(1).FaceColor = standard_colors('pbo');
b(2).FaceColor = standard_colors('etc');

errorbar(b(1).XEndPoints, YData(:,1), EData(:,1), 'k','LineStyle','none','CapSize',2);
errorbar(b(2).XEndPoints, YData(:,2), EData(:,2), 'k','LineStyle','none','CapSize',2);

text(Ax(ai), 1, 56, '*', 'FontSize', 8, 'HorizontalAlignment', 'center');

l = legend(Ax(ai), b, {'Placebo','THC/CBD'}, 'Box','off','FontSize',7);
l.Position = [0.28 0.8125 0.15 0.1];
l.ItemTokenSize = [5 5];

Ax(ai).FontSize = 8;
Ax(ai).TickLength = [0 0];
Ax(ai).XTick = 1:3;
Ax(ai).XTickLabel = {'N1','N2','N3'};
Ax(ai).YTick = 0:10:100;
Ax(ai).YLabel.String = 'Arousal Index (N/hr)';
Ax(ai).YLabel.FontSize = 8;
Ax(ai).OuterPosition = [0 0.5 0.5 0.46];
Ax(ai).Color = [0.95 0.96 0.98];
Ax(ai).Box = 'on';
Ax(ai).YGrid = 'on';

Ax(ai).Title.String = 'A';
Ax(ai).Title.FontSize = 8;
Ax(ai).Title.HorizontalAlignment = 'left';
Ax(ai).Title.Units = 'normalized';
Ax(ai).Title.Position(1:2) = [-0.22 1.1];

% Panel B — Histogram of Arousal Durations
ai = ai + 1;
Ax(ai) = axes('FontSize',8);

idx = ARO.stage <= -1 & strcmpi(ARO.aro_type, '1arousal');

YData = [ ...
    histcounts(ARO.duration(idx & ARO.cond=="placebo"), BinEdges)', ...
    histcounts(ARO.duration(idx & ARO.cond=="etc120"), BinEdges)'];

b = bar(Ax(ai), BinCenters, YData, 1, 'GroupWidth',0.9);
b(1).FaceColor = standard_colors('pbo');
b(2).FaceColor = standard_colors('etc');

Ax(ai).FontSize = 8;
Ax(ai).TickLength = [0 0];
Ax(ai).Color = [0.95 0.96 0.98];
Ax(ai).Box = 'on';
Ax(ai).XGrid = 'on';
Ax(ai).XTick = [0 round(BinCenters(2:2:end))];
Ax(ai).XTickLabelRotation = 0;
Ax(ai).XLim = [0 32];
Ax(ai).YTick = [0 max(YData(:))];
Ax(ai).XLabel.String = 'Duration (s)';
Ax(ai).XLabel.FontSize = 8;
Ax(ai).YLabel.String = 'Frequency (N)';
Ax(ai).YLabel.FontSize = 8;
Ax(ai).YLabel.Position(1) = -6;

Ax(ai).Title.String = 'B';
Ax(ai).Title.FontSize = 8;
Ax(ai).Title.HorizontalAlignment = 'left';
Ax(ai).Title.Units = 'normalized';
Ax(ai).Title.Position(1) = -0.22;

Ax(ai).Position = Ax(ai-1).Position + [0 -0.45 0 -0.05];

% Panels C–E — Probability of Awakening by Duration for N1, N2, N3
stages = [-1, -2, -3];
labels = {'C   NREM stage-1', 'D   NREM stage-2', 'E   NREM stage-3'};
positions = [0.63, 0.31, 0.00];

for s = 1:length(stages)
    ai = ai + 1;
    Ax(ai) = axes('NextPlot','add');

    idx = ARO.stage==stages(s) & strcmpi(ARO.aro_type,'1arousal');
    ARO.duration_bin = discretize(ARO.duration, BinEdges);

    P = groupsummary(ARO(idx,:), {'participant_id','cond'}, 'all', {'is_awakening','duration'}); %#ok<NASGU>

    G = groupsummary(ARO(idx,:), {'duration_bin','cond'}, 'all', {'is_awakening','duration'});
    G.PairedGroupCount = [G.GroupCount(1:end-1) + G.GroupCount(2:end); nan];
    G.PairedGroupCount(2:2:end) = G.PairedGroupCount(1:2:end-1);

    G.conf_int_duration = tinv(0.975, G.PairedGroupCount-1).*G.std_duration./sqrt(G.PairedGroupCount);
    G.conf_is_awakening = tinv(0.975, G.PairedGroupCount-1).*G.std_is_awakening./sqrt(G.PairedGroupCount);

    plot([3.5, 3.5], [0, 1.05], '-k')

    errorbar(1.1, S(s).pred(1), S(s).pred(1) - S(s).ci(1, 1), 'o', 'MarkerSize', 3, 'Color', standard_colors('black'), 'MarkerEdgeColor', standard_colors('black'), 'MarkerFaceColor', standard_colors('pbo'), 'CapSize', 1.5)
    errorbar(2.4, S(s).pred(2), S(s).pred(2) - S(s).ci(2, 1), 'o', 'MarkerSize', 3, 'Color', standard_colors('black'), 'MarkerEdgeColor', standard_colors('black'), 'MarkerFaceColor', standard_colors('etc'), 'CapSize', 1.5)

    errorpatch(Ax(ai), G.mean_duration(G.cond=="placebo"), G.mean_is_awakening(G.cond=="placebo"), G.conf_is_awakening(G.cond=="placebo"), ...
        'FaceAlpha', 0.33, ...
        'FaceColor', standard_colors('pbo'));

    errorpatch(Ax(ai), G.mean_duration(G.cond=="etc120"), G.mean_is_awakening(G.cond=="etc120"), G.conf_is_awakening(G.cond=="etc120"), ...
        'FaceAlpha', 0.33, ...
        'FaceColor', standard_colors('etc'));

    % Plot
    plot(G.mean_duration(G.cond=="placebo"), ...
        G.mean_is_awakening(G.cond=="placebo"), ...
        '-', 'Color', standard_colors('pbo'), 'LineWidth', 1.5);

    plot(G.mean_duration(G.cond=="etc120"), ...
        G.mean_is_awakening(G.cond=="etc120"), ...
        '-', 'Color', standard_colors('etc'), 'LineWidth', 1.5);

    Ax(ai).TickLength = [0 0];
    Ax(ai).FontSize = 8;
    Ax(ai).Color = [0.95 0.96 0.98];
    Ax(ai).Box = 'on';

    Ax(ai).YLim = [0 1.05];
    Ax(ai).YTick = [0 1];
    Ax(ai).YLabel.String = 'Pr(Awake|Ar)';
    Ax(ai).YLabel.FontSize = 8;

    Ax(ai).XLim = [0 32];
    Ax(ai).XTick = [round(BinCenters(2:2:end))];
    Ax(ai).XGrid = 'on';

    if s==3
        Ax(ai).XLabel.String = 'Duration (s)';
        Ax(ai).XLabel.FontSize = 8;
    end

    Ax(ai).Title.String = ['\bf', labels{s}(1), '\rm           ', labels{s}(2:end)];
    Ax(ai).Title.FontSize = 8;
    Ax(ai).Title.HorizontalAlignment = 'left';
    Ax(ai).Title.Units = 'normalized';
    Ax(ai).Title.Position(1) = -0.22;

    Ax(ai).Position = [0.5 positions(s) 0.5 0.33] + [0.1 0.11 -0.12 -0.12];
end

exportgraphics(Fig, 'figures/fig_arousal_descriptives.png', 'Resolution', 600);

end
