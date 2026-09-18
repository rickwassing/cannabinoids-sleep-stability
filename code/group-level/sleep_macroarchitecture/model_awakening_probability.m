function [ARO, S] = model_awakening_probability(ARO)
% -------------------------------------------------------------------------
% Mixed-effects modelling: Onset, Duration, Awakening probability
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
stats = compare(mdl_null, mdl_base) %#ok<NOPRT,NASGU>
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

% -------------------------------------------------------------------------
% Stage-specific binomial models
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
    [S(i).pred, S(i).ci] = calc_risk_ratio(ARO, S(i).mdl);
    S(i).T.Factor = {'THC/CBD'; 'Arousal duration'};
    S(i).T.Reference = {'Placebo'; ''};
    S(i).T.Beta = arrayfun(@(v) sprintf('%.2f', v), S(i).mdl.Coefficients.Estimate(2:end), 'UniformOutput', false);
    S(i).T.SE = arrayfun(@(v) sprintf('%.2f', v), S(i).mdl.Coefficients.SE(2:end), 'UniformOutput', false);
    S(i).T.tStat = arrayfun(@(v) sprintf('%.2f', v), S(i).mdl.Coefficients.tStat(2:end), 'UniformOutput', false);
    S(i).T.pValue = arrayfun(@(v) stringify_pvalue(v), S(i).mdl.Coefficients.pValue(2:end), 'UniformOutput', false);
    S(i).T.Exp_Beta = arrayfun(@(v) sprintf('%.2f', v), exp(S(i).mdl.Coefficients.Estimate(2:end)), 'UniformOutput', false);
    S(i).T.CI = arrayfun(@(l, h) sprintf('[%.2f - %.2f]', l, h), exp(S(i).mdl.Coefficients.Lower(2:end)), exp(S(i).mdl.Coefficients.Upper(2:end)), 'UniformOutput', false);

    S(i).coefTable %#ok<NOPRT>

end

end
