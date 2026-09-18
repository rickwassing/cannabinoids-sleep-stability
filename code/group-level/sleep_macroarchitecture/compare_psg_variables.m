function T = compare_psg_variables(PSG)
% -------------------------------------------------------------------------
% Compare PSG variables (Table 1) between conditions using linear
% mixed-effects models, and write the descriptive/inferential summary to
% 'figures/tbl_psg_descriptives.csv'.
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
    T.pval{i} = stringify_pvalue(mdl.Coefficients.pValue(2));
end

writetable(T, 'figures/tbl_psg_descriptives.csv');

end
