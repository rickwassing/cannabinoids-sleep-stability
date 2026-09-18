PSG = readtable('./phenotype/2024-07-26T171508_psg-variables.csv');

Fig = figure();
Fig.Units = 'centimeters';
Fig.Position = [10 10 8.5 8.5];

Ax = axes('NextPlot', 'add');
Ax.Color = [0.96 0.97 0.99];
Ax.Box = 'on';
Ax.YGrid = 'on';
Ax.FontSize = 8;
Ax.TickLength = [0, 0];
Ax.XLim = [0.33, 2.67];
Ax.XTick = [];
Ax.YLim = [0 480];
Ax.YTick = 0:60:480;
Ax.YLabel.String = 'Total sleep time (min)';
Ax.YLabel.FontSize = 10;
Ax.Title.String = '–24 min, p < 0.05';
Ax.Title.FontSize = 8;
Ax.Title.FontWeight = 'normal';

XData = 1:2;
BData = [...
    mean(PSG.tst(strcmpi(PSG.condition, 'placebo'))), 0; ...
    0, mean(PSG.tst(strcmpi(PSG.condition, 'etc120')))];
YData = [...
    mean(PSG.tst(strcmpi(PSG.condition, 'placebo'))), ...
    mean(PSG.tst(strcmpi(PSG.condition, 'etc120')))];
EData = [...
    std(PSG.tst(strcmpi(PSG.condition, 'placebo'))), ...
    std(PSG.tst(strcmpi(PSG.condition, 'etc120')))];
h = bar(Ax, BData, 'stacked');
errorbar(Ax, XData, YData, EData, 'k', 'LineStyle', 'none')

h(1).FaceColor = standard_colors('green');
h(2).FaceColor = standard_colors('orange');

Ax.OuterPosition = [0 0.5 0.4 0.5];

l = legend(h, {'PBO', 'THC/CBD'});
l.Box = 'off';
l.ItemTokenSize = [10, 10];
l.Position = [0.85 Ax.Position(2)+Ax.Position(4)-0.1 0.1 0.1];
l.String = {'PBO', 'THC/CBD'};

% -------------------------------------------------------------------------

Ax = axes('NextPlot', 'add');
Ax.Color = [0.96 0.97 0.99];
Ax.Box = 'on';
Ax.YGrid = 'on';
Ax.FontSize = 8;
Ax.TickLength = [0, 0];
Ax.XLim = [0.33, 2.67];
Ax.XTick = [];
Ax.YLim = [0 60];
Ax.YTick = 0:10:100;
Ax.YLabel.String = 'NREM stage 2 (%)';
Ax.YLabel.FontSize = 10;
Ax.Title.String = '+5%, p < 0.02';
Ax.Title.FontSize = 8;
Ax.Title.FontWeight = 'normal';

XData = 1:2;
YData = [...
    mean(PSG.n2_pct(strcmpi(PSG.condition, 'placebo'))), ...
    mean(PSG.n2_pct(strcmpi(PSG.condition, 'etc120')))];
EData = [...
    std(PSG.n2_pct(strcmpi(PSG.condition, 'placebo'))), ...
    std(PSG.n2_pct(strcmpi(PSG.condition, 'etc120')))];
h = bar(Ax, XData, YData);
errorbar(Ax, XData, YData, EData, 'k', 'LineStyle', 'none')

h.FaceColor = 'flat';
h.CData(1, :) = standard_colors('green');
h.CData(2, :) = standard_colors('orange');

Ax.OuterPosition = [0.4 0.5 0.4 0.5];

% -------------------------------------------------------------------------

Ax = axes('NextPlot', 'add');
Ax.Color = [0.96 0.97 0.99];
Ax.Box = 'on';
Ax.YGrid = 'on';
Ax.FontSize = 8;
Ax.TickLength = [0, 0];
Ax.XLim = [0.33, 2.67];
Ax.XTick = [];
Ax.YLim = [0 30];
Ax.YTick = 0:5:300;
Ax.YLabel.String = 'REM (%)';
Ax.YLabel.FontSize = 10;
Ax.Title.String = '–8%, p < 0.01';
Ax.Title.FontSize = 8;
Ax.Title.FontWeight = 'normal';

XData = 1:2;
YData = [...
    mean(PSG.rem_pct(strcmpi(PSG.condition, 'placebo'))), ...
    mean(PSG.rem_pct(strcmpi(PSG.condition, 'etc120')))];
EData = [...
    std(PSG.rem_pct(strcmpi(PSG.condition, 'placebo'))), ...
    std(PSG.rem_pct(strcmpi(PSG.condition, 'etc120')))];
h = bar(Ax, XData, YData);
errorbar(Ax, XData, YData, EData, 'k', 'LineStyle', 'none')

h.FaceColor = 'flat';
h.CData(1, :) = standard_colors('green');
h.CData(2, :) = standard_colors('orange');

Ax.OuterPosition = [0 0 0.4 0.5];

% -------------------------------------------------------------------------

Ax = axes('NextPlot', 'add');
Ax.Color = [0.96 0.97 0.99];
Ax.Box = 'on';
Ax.YGrid = 'on';
Ax.FontSize = 8;
Ax.TickLength = [0, 0];
Ax.XLim = [0.33, 2.67];
Ax.XTick = [];
Ax.YLim = [0 270];
Ax.YTick = 0:30:300;
Ax.YLabel.String = 'REML (min)';
Ax.YLabel.FontSize = 10;
Ax.Title.String = '+60 min, p < 0.01';
Ax.Title.FontSize = 8;
Ax.Title.FontWeight = 'normal';

XData = 1:2;
YData = [...
    mean(PSG.remlat(strcmpi(PSG.condition, 'placebo'))), ...
    mean(PSG.remlat(strcmpi(PSG.condition, 'etc120')))];
EData = [...
    std(PSG.remlat(strcmpi(PSG.condition, 'placebo'))), ...
    std(PSG.remlat(strcmpi(PSG.condition, 'etc120')))];
h = bar(Ax, XData, YData);
errorbar(Ax, XData, YData, EData, 'k', 'LineStyle', 'none')

h.FaceColor = 'flat';
h.CData(1, :) = standard_colors('green');
h.CData(2, :) = standard_colors('orange');

Ax.OuterPosition = [0.4 0 0.4 0.5];

exportgraphics(Fig, './figures/fig_psgbargraphs.png', 'Resolution', 600)

%%

mdl = fitlme(PSG, 'i_aro_rem ~ 1 + condition + (1|participant_id)')