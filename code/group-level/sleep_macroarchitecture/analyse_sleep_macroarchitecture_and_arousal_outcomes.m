function analyse_sleep_macroarchitecture_and_arousal_outcomes()
% Load data
PSG = readtable('phenotype/2024-07-26T171508_psg-variables.csv');
ARO = readtable('phenotype/2024-07-26T171508_arousals.csv');
ARO.duration(ARO.duration > 30) = 30;  % Cap extreme durations

% Stage coding
ARO.stage_str = cell(size(ARO,1),1);
ARO.stage_str(ARO.stage ==  1) = {'r'};
ARO.stage_str(ARO.stage ==  0) = {'w'};
ARO.stage_str(ARO.stage == -1) = {'n1'};
ARO.stage_str(ARO.stage == -2) = {'n2'};
ARO.stage_str(ARO.stage == -3) = {'n3'};

% ------------------------------------------------------------------------
%  Compare PSG variables between conditions
% -------------------------------------------------------------------------
vars = { ...
    'tst','se','sol','waso','n1_pct','n2_pct','n3_pct','rem_pct','remlat', ...
    'n_aro_nrem','n_aro_n1','n_aro_n2','n_aro_n3','n_aro_rem', ...
    'i_aro_nrem','i_aro_n1','i_aro_n2','i_aro_n3','i_aro_rem'};

idxETC = strcmpi(PSG.condition, 'etc120');

% Preallocate descriptive table
T = table();
T.Variable   = vars';
T.EtcMeanSD  = cell(length(vars),1);
T.EtcRange   = cell(length(vars),1);
T.PboMeanSD  = cell(length(vars),1);
T.PboRange   = cell(length(vars),1);
T.tstat      = cell(length(vars),1);
T.pval       = cell(length(vars),1);

for i = 1:length(vars)
    v = PSG.(vars{i});

    % Descriptives
    T.EtcMeanSD{i} = sprintf('%.1f (%.1f)', mean(v(idxETC)), std(v(idxETC)));
    T.PboMeanSD{i} = sprintf('%.1f (%.1f)', mean(v(~idxETC)), std(v(~idxETC)));

    T.EtcRange{i} = sprintf('[%.1f - %.1f]', min(v(idxETC)), max(v(idxETC)));
    T.PboRange{i} = sprintf('[%.1f - %.1f]', min(v(~idxETC)), max(v(~idxETC)));

    % Linear mixed model
    mdl = fitlme(PSG, sprintf('%s ~ 1 + condition + (1|participant_id)', vars{i}), ...
        'DummyVarCoding', 'effects');
    T.tstat{i} = sprintf('%.2f', mdl.Coefficients.tStat(2));
    T.DFE{i} = sprintf('%.i', mdl.DFE);
    T.pval{i} = stringifypvalue(mdl.Coefficients.pValue(2));
end

writetable(T, 'figures/tbl_psg_descriptives.csv');

%% ------------------------------------------------------------------------
%  Mixed-effects modelling: Onset, Duration, Awakening probability
% -------------------------------------------------------------------------
clc

idxAll = ARO.stage <= -1 & strcmpi(ARO.aro_type,'1arousal');

% Z-Score duration for statistical modeling
ARO.duration_z = zscore(ARO.duration);
ARO.duration_z(idxAll) = zscore(ARO.duration(idxAll));

% Models
disp('#################################################################')
disp('DURATION MODEL')
mdl_dur = fitglme(ARO(idxAll,:), 'duration ~ 1 + cond*stage_str + (1|stage_str) + (1|participant_id)', ...
    'Distribution','Normal','DummyVarCoding','effects');
disp(mdl_dur);
disp(anova(mdl_dur));

disp('#################################################################')
disp('AWAKENING MODEL')
disp('-----------------------------------------------------------------')
disp('Null models')
mdl_null = fitglme(ARO(idxAll,:), 'is_awakening ~ 1', ...
    'Distribution','Binomial','DummyVarCoding','effects', 'FitMethod','Laplace');
disp(mdl_null);
disp(anova(mdl_null));
disp('-----------------------------------------------------------------')
disp('Base model')
mdl_base = fitglme(ARO(idxAll,:), 'is_awakening ~ 1 + (1|participant_id) + (1|stage_str)', ...
    'Distribution','Binomial','DummyVarCoding','effects', 'FitMethod','Laplace');
disp(mdl_base);
disp(anova(mdl_base));
stats = compare(mdl_null, mdl_base)
disp('-----------------------------------------------------------------')
disp('Full model')
mdl_full = fitglme(ARO(idxAll,:), 'is_awakening ~ 1 + cond*stage_str*duration_z + (1|stage_str) + (1|participant_id)', ...
    'Distribution','Binomial','DummyVarCoding','effects', 'FitMethod','Laplace');
disp(mdl_full);
disp(anova(mdl_full));
disp('-----------------------------------------------------------------')
disp('Adapted model')
mdl_fin = fitglme(ARO(idxAll,:), 'is_awakening ~ 1 + cond + stage_str + duration_z + stage_str:duration_z + cond:stage_str + (1|stage_str) + (1|participant_id)', ...
    'Distribution','Binomial','DummyVarCoding','effects', 'FitMethod','Laplace');
coefTable = mdl_fin.Coefficients;
coefTable.OddsRatio = exp(coefTable.Estimate);
coefTable.OddsLo = exp(coefTable.Lower);
coefTable.OddsHi = exp(coefTable.Upper);
disp(coefTable);
disp(anova(mdl_fin));

%% ------------------------------------------------------------------------
%  Stage-specific binomial models
% -------------------------------------------------------------------------
clc

T = table();

S = struct();
S.formula = '';
S.mdl = [];
S.pred = [];
S.ci = [];
S.T = T;
S = repmat(S, 3, 1);

stages = [-1, -2, -3];
for i = 1:length(stages)

    idx = ARO.stage == stages(i) & strcmpi(ARO.aro_type,'1arousal');
    ARO.duration_z(idx) = zscore(ARO.duration(idx));
    S(i).mdl = fitglme(ARO(idx, :), 'is_awakening ~ 1 + cond + duration_z + (1|participant_id)', ...
        'Distribution','Binomial');
    S(i).coefTable = S(i).mdl.Coefficients;
    S(i).coefTable.OddsRatio = exp(S(i).coefTable.Estimate);
    S(i).coefTable.OddsLo = exp(S(i).coefTable.Lower);
    S(i).coefTable.OddsHi = exp(S(i).coefTable.Upper);
    S(i).formula = char(S(i).mdl.Formula);
    [S(i).pred, S(i).ci] = calcRiskRatio(ARO, S(i).mdl);
    S(i).T.Factor = {'THC/CBD'; 'Arousal duration'};
    S(i).T.Reference = {'Placebo'; ''};
    S(i).T.Beta = arrayfun(@(v) sprintf('%.2f', v), S(i).mdl.Coefficients.Estimate(2:end), 'UniformOutput', false);
    S(i).T.SE = arrayfun(@(v) sprintf('%.2f', v), S(i).mdl.Coefficients.SE(2:end), 'UniformOutput', false);
    S(i).T.tStat = arrayfun(@(v) sprintf('%.2f', v), S(i).mdl.Coefficients.tStat(2:end), 'UniformOutput', false);
    S(i).T.pValue = arrayfun(@(v) stringifypvalue(v), S(i).mdl.Coefficients.pValue(2:end), 'UniformOutput', false);
    S(i).T.Exp_Beta = arrayfun(@(v) sprintf('%.2f', v), exp(S(i).mdl.Coefficients.Estimate(2:end)), 'UniformOutput', false);
    S(i).T.CI = arrayfun(@(l, h) sprintf('[%.2f - %.2f]', l, h), exp(S(i).mdl.Coefficients.Lower(2:end)), exp(S(i).mdl.Coefficients.Upper(2:end)), 'UniformOutput', false);

    S(i).coefTable

end

%% ------------------------------------------------------------------------
%  Figure: Arousal descriptives
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
b(1).FaceColor = css_standard_colors('pbo');
b(2).FaceColor = css_standard_colors('etc');

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
b(1).FaceColor = css_standard_colors('pbo');
b(2).FaceColor = css_standard_colors('etc');

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

    P = groupsummary(ARO(idx,:), {'participant_id','cond'}, 'all', {'is_awakening','duration'});

    G = groupsummary(ARO(idx,:), {'duration_bin','cond'}, 'all', {'is_awakening','duration'});
    G.PairedGroupCount = [G.GroupCount(1:end-1) + G.GroupCount(2:end); nan];
    G.PairedGroupCount(2:2:end) = G.PairedGroupCount(1:2:end-1);

    G.conf_int_duration = tinv(0.975, G.PairedGroupCount-1).*G.std_duration./sqrt(G.PairedGroupCount);
    G.conf_is_awakening = tinv(0.975, G.PairedGroupCount-1).*G.std_is_awakening./sqrt(G.PairedGroupCount);

    plot([3.5, 3.5], [0, 1.05], '-k')

    errorbar(1.1, S(s).pred(1), S(s).pred(1) - S(s).ci(1, 1), 'o', 'MarkerSize', 3, 'Color', css_standard_colors('black'), 'MarkerEdgeColor', css_standard_colors('black'), 'MarkerFaceColor', css_standard_colors('pbo'), 'CapSize', 1.5)
    errorbar(2.4, S(s).pred(2), S(s).pred(2) - S(s).ci(2, 1), 'o', 'MarkerSize', 3, 'Color', css_standard_colors('black'), 'MarkerEdgeColor', css_standard_colors('black'), 'MarkerFaceColor', css_standard_colors('etc'), 'CapSize', 1.5)

    errorpatch(Ax(ai), G.mean_duration(G.cond=="placebo"), G.mean_is_awakening(G.cond=="placebo"), G.conf_is_awakening(G.cond=="placebo"), ...
        'FaceAlpha', 0.33, ...
        'FaceColor', css_standard_colors('pbo'));

    errorpatch(Ax(ai), G.mean_duration(G.cond=="etc120"), G.mean_is_awakening(G.cond=="etc120"), G.conf_is_awakening(G.cond=="etc120"), ...
        'FaceAlpha', 0.33, ...
        'FaceColor', css_standard_colors('etc'));

    % Plot
    plot(G.mean_duration(G.cond=="placebo"), ...
        G.mean_is_awakening(G.cond=="placebo"), ...
        '-', 'Color', css_standard_colors('pbo'), 'LineWidth', 1.5);

    plot(G.mean_duration(G.cond=="etc120"), ...
        G.mean_is_awakening(G.cond=="etc120"), ...
        '-', 'Color', css_standard_colors('etc'), 'LineWidth', 1.5);

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


function str = stringifypvalue(p)
if p < 1e-4
    str = sprintf('<10^%.0f', round(log10(p) + 1));
elseif p < 0.001
    str = sprintf('%.5f', p);
elseif p < 0.01
    str = sprintf('%.4f', p);
else
    str = sprintf('%.3f', p);
end
end

function [predProb, predCI] = calcRiskRatio(ARO, mdl)
% Predictions
predTbl = ARO(1:2,:);
predTbl.cond(1) = {'placebo'};
predTbl.cond(2) = {'etc120'};
predTbl.duration = [0; 0];
predTbl.duration_z = [0; 0];

[predProb, predCI] = predict(mdl, predTbl, 'Conditional', false);

end