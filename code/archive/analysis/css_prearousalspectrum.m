clear
clc

[~, this_host] = system('hostname');
this_host = strtrim(this_host);

switch this_host
    case 'wimr-HP-Z6'
        stage = {'aroboutrem'};
    case 'WIMR-HPZ8-01'
        stage = {'aroboutnrem2'};
    case 'Ricks-MacBook-Air.local'
        stage = {'aroboutnrem2'};
    case 'mnc.local'
        stage = {'aroboutrem'};
    otherwise
        stage = {'aroboutnrem2'};
end

Files = dir(sprintf('derivatives/aroepochs/sub-*/ses-*/sub*%s*.mat', stage{:}));

AData = struct();
AData.cs = nan(0, 241);
AData.aw = nan(0, 241);

CData = struct();
CData.cs = nan(178, 241, 30, 0, 'single');
CData.aw = nan(178, 241, 30, 0, 'single');

Subj = struct();
Subj.cs = [];
Subj.aw = [];

Cond = struct();
Cond.cs = [];
Cond.aw = [];

smooth_factor = 6; % seconds
do_smooth = true;
do_normalize = true;
rt = now();
for i = 1:length(Files)
    kv = filename2struct(Files(i).name);
    % Load
    EEG = load(fullfile(Files(i).folder, Files(i).name));
    ARO = EEG.ARO;
    clear EEG;
    % Normalize
    if do_normalize
        idx_t = ARO.times < 0; %#ok<UNRCH>
        freqstep = mean(diff(ARO.specfreqs));
        grandtotalpower = squeeze(mean(mean(mean(sum(ARO.data(:, idx_t, :, :), 3, 'omitnan') .* freqstep, 1, 'omitnan'), 2, 'omitnan'), 4, 'omitnan'));
        ARO.data = (ARO.data.*freqstep) ./ grandtotalpower;
    end
    % Average across frequencies (whole freqs)
    d = nan(size(ARO.data, 1), length(ARO.times), 30, size(ARO.data, 4), 'single');
    for f = 1:max(ARO.specfreqs)+0.001
        idx_f = find(ARO.specfreqs > (f-1+0.001) & ARO.specfreqs <= f+0.001);
        d(:, :, f, :) = mean(ARO.data(:, :, idx_f, :), 3, 'omitnan');
    end
    ARO.data = d;
    if do_smooth
        % Smooth time
        srate = 1/mean(diff(ARO.times));
        for c = 1:size(ARO.data, 1)
            for f = 1:size(ARO.data, 3)
                for e = 1:size(ARO.data, 4)
                    ARO.data(c, :, f, e) = smooth(ARO.data(c, :, f, e), smooth_factor * srate + 1);
                end
            end
        end
    end
    % -------------------------------------------------------------------------
    % Extract awakening and continued sleep trials
    idx_aw = ARO.event.is_awakening;
    AData.cs = cat(1, AData.cs, ARO.aro(~idx_aw, :));
    AData.aw = cat(1, AData.aw, ARO.aro(idx_aw, :));
    % -------------------------------------------------------------------------
    CData.cs = cat(4, CData.cs, ARO.data(:, :, :, ~idx_aw));
    CData.aw = cat(4, CData.aw, ARO.data(:, :, :, idx_aw));
    % -------------------------------------------------------------------------
    Subj.cs = [Subj.cs; repmat({kv.sub}, sum(~idx_aw), 1)];
    Subj.aw = [Subj.aw; repmat({kv.sub}, sum(idx_aw), 1)];
    % -------------------------------------------------------------------------
    Cond.cs = [Cond.cs; repmat({kv.ses}, sum(~idx_aw), 1)];
    Cond.aw = [Cond.aw; repmat({kv.ses}, sum(idx_aw), 1)];
    % -------------------------------------------------------------------------
    % Remaining time
    rt = remainingTime(rt, length(Files));
end
% -------------------------------------------------------------------------
% Save
foutname = sprintf('group-level/%s_prearo-%s.mat', datestr(now, 'yyyymmddTHHMM'), stage{:});
save(foutname, 'AData', '-mat', '-v7.3');
save(foutname, 'CData', '-append');
save(foutname, 'Cond', '-append');
save(foutname, 'Subj', '-append');
save(foutname, 'do_normalize', '-append');
save(foutname, 'do_smooth', '-append');
save(foutname, 'smooth_factor', '-append');

%% Test for differences in N2 sigma power between ETC120 and placebo
idx_cond_cs = strcmpi(Cond.cs, 'etc120');
idx_cond_aw = strcmpi(Cond.aw, 'etc120');

specfreqs = 1:30;
times = -90:0.5:30;
idx_t = times < -4;
idx_f = specfreqs >= 10.9 & specfreqs <= 16.1;

AvCDataPLC = squeeze(mean(mean(cat(4, CData.cs(:, idx_t, idx_f, ~idx_cond_cs), CData.aw(:, idx_t, idx_f, ~idx_cond_aw)), 2, 'omitnan'), 3, 'omitnan'));
AvCDataETC = squeeze(mean(mean(cat(4, CData.cs(:, idx_t, idx_f, idx_cond_cs), CData.aw(:, idx_t, idx_f, idx_cond_aw)), 2, 'omitnan'), 3, 'omitnan'));

tbl = table();
tbl.y = nan(size(AvCDataPLC, 2) + size(AvCDataETC, 2), 1);
tbl.cond = [-1.*ones(size(AvCDataPLC, 2), 1); 1.*ones(size(AvCDataETC, 2), 1)];
tbl.subj = [Subj.cs(~idx_cond_cs); Subj.aw(~idx_cond_aw); Subj.cs(idx_cond_cs); Subj.aw(idx_cond_aw)];

tstat = nan(178, 1);
for i = 1:178
    tbl.y = double([AvCDataPLC(i, :), AvCDataETC(i, :)]');
    tbl.y = zscore(tbl.y);
    mdl = fitlme(tbl, 'y ~ 1 + cond + (1|subj)');
    tstat(i) = mdl.Coefficients.tStat(2);
end

%% APPLY GLM
% -------------------------------------------------------------------------
% Get channel location clusters
ARO.chanlocs = channel_clusters(ARO.chanlocs, 'MFF');
% -------------------------------------------------------------------------
% Init table for modelling
tbl = table();
tbl.is_awake = [true(size(CData.aw, 4), 1); false(size(CData.cs, 4), 1)];
tbl.cond = [Cond.aw; Cond.cs];
tbl.sub = [Subj.aw; Subj.cs];
% -------------------------------------------------------------------------
% Output structure
GLM = struct();
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Linear mixed model, POW ~ 1 + IS_AWAKE*COND
GLM.logpow.awake.tstat = nan(6, 180, 30);
GLM.logpow.awake.pval = nan(6, 180, 30);
GLM.logpow.cond.tstat = nan(6, 180, 30);
GLM.logpow.cond.pval = nan(6, 180, 30);
GLM.logpow.intx.tstat = nan(6, 180, 30);
GLM.logpow.intx.pval = nan(6, 180, 30);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Linear mixed model, POW ~ 1 + IS_AWAKE + COND
GLM.rem.awake.tstat = nan(6, 180, 30);
GLM.rem.awake.pval = nan(6, 180, 30);
GLM.rem.cond.tstat = nan(6, 180, 30);
GLM.rem.cond.pval = nan(6, 180, 30);
% -------------------------------------------------------------------------
% Output structure
PostHoc = struct();
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Linear mixed model, POW ~ 1 + IS_AWAKE (separate for placebo and ETC120)
PostHoc.plc.awake.tstat = nan(6, 180, 30);
PostHoc.plc.awake.pval = nan(6, 180, 30);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
PostHoc.etc.awake.tstat = nan(6, 180, 30);
PostHoc.etc.awake.pval = nan(6, 180, 30);

% -------------------------------------------------------------------------
rt = now(); %#ok<TNOW1>
for c = 1:6
    % ---------------------------------------------------------------------
    % Get channel indices for this cluster
    idx_c = [ARO.chanlocs.cluster] == c;
    for t = 1:180
        for f = 1:30
            % -------------------------------------------------------------
            % Replace power values in the table
            tbl.pow = double([...
                squeeze(mean(mean(mean(CData.aw(idx_c, t, f, :), 3, 'omitnan'), 2, 'omitnan'), 1, 'omitnan')); ...
                squeeze(mean(mean(mean(CData.cs(idx_c, t, f, :), 3, 'omitnan'), 2, 'omitnan'), 1, 'omitnan'))]);
            tmp_aw = squeeze(mean(mean(mean(CData.aw(idx_c, t, f, :), 3, 'omitnan'), 2, 'omitnan'), 1, 'omitnan'));
            tmp_cs = squeeze(mean(mean(mean(CData.cs(idx_c, t, f, :), 3, 'omitnan'), 2, 'omitnan'), 1, 'omitnan'));
            tmp_aw(tmp_aw <= 0) = min([tmp_aw(tmp_aw > 0); tmp_cs(tmp_cs > 0)]);
            tmp_cs(tmp_cs <= 0) = min([tmp_aw(tmp_aw > 0); tmp_cs(tmp_cs > 0)]);
            tbl.log_mean_pow = double([...
                log10(tmp_aw); ...
                log10(tmp_cs)]);
            % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            % Zscore so the maths behave normally
            tbl.pow = nanzscore(tbl.pow);
            tbl.log_mean_pow = nanzscore(tbl.log_mean_pow);
            % -------------------------------------------------------------
            % Fit full linear mixed model
            glme = fitlme(tbl, 'log_mean_pow ~ 1 + is_awake*cond + (1|sub)', 'DummyVarCoding', 'effects');
            % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            GLM.logpow.awake.tstat(c, t, f) = glme.Coefficients.tStat(2);
            GLM.logpow.awake.pval(c, t, f) = glme.Coefficients.pValue(2);
            GLM.logpow.cond.tstat(c, t, f) = glme.Coefficients.tStat(3);
            GLM.logpow.cond.pval(c, t, f) = glme.Coefficients.pValue(3);
            GLM.logpow.intx.tstat(c, t, f) = glme.Coefficients.tStat(4);
            GLM.logpow.intx.pval(c, t, f) = glme.Coefficients.pValue(4);
            % -------------------------------------------------------------
            % Fit simple linear mixed model
            glme_main = fitlme(tbl, 'log_mean_pow ~ 1 + is_awake + cond + (1|sub)', 'DummyVarCoding', 'effects');
            % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            GLM.rem.awake.tstat(c, t, f) = glme_main.Coefficients.tStat(2);
            GLM.rem.awake.pval(c, t, f) = glme_main.Coefficients.pValue(2);
            GLM.rem.cond.tstat(c, t, f) = glme_main.Coefficients.tStat(3);
            GLM.rem.cond.pval(c, t, f) = glme_main.Coefficients.pValue(3);
            % -------------------------------------------------------------
            % Fit linear mixed model for placebo
            idx_rows = strcmpi(tbl.cond, 'placebo');
            glme_plc = fitlme(tbl(idx_rows, :), 'log_mean_pow ~ 1 + is_awake + (1|sub)', 'DummyVarCoding', 'effects');
            % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            PostHoc.plc.awake.tstat(c, t, f) = glme_plc.Coefficients.tStat(2);
            PostHoc.plc.awake.pval(c, t, f) = glme_plc.Coefficients.pValue(2);
            % -------------------------------------------------------------
            % Fit linear mixed model
            idx_rows = strcmpi(tbl.cond, 'etc120');
            glme_etc = fitlme(tbl(idx_rows, :), 'log_mean_pow ~ 1 + is_awake + (1|sub)', 'DummyVarCoding', 'effects');
            % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            PostHoc.etc.awake.tstat(c, t, f) = glme_etc.Coefficients.tStat(2);
            PostHoc.etc.awake.pval(c, t, f) = glme_etc.Coefficients.pValue(2);
            % -------------------------------------------------------------
            % Display remaining time
            rt = remainingTime(rt, 6*180*30, 'simple');
        end
    end
end

% RUN THIS AGAIN TO SAVE 'POSTHOC'
save(sprintf('group-level/%s_modeldata-%s.mat', datestr(now, 'yyyymmddTHHMM'), stage{:}), 'GLM', 'PostHoc', 'glme', 'glme_main', 'glme_plc', 'glme_etc', '-mat', '-v7.3') %#ok<DATST,TNOW1>

%% RUN POST-HOC ANALYSIS USING PALM
% -------------------------------------------------------------------------
% NREM SLEEP FINDINGS
% -------------------------------------------------------------------------
% Create dirs if needed
if exist('group-level/lme_prearousalsigma_placebo', 'dir') == 0
    mkdir('group-level/lme_prearousalsigma_placebo/input')
    mkdir('group-level/lme_prearousalsigma_placebo/output')
    copyfile('group-level/paired-t-test_clustermass_perm/input/surface.srf', 'group-level/lme_prearousalsigma_placebo/input/')
    copyfile('group-level/paired-t-test_clustermass_perm/input/chanlocs.sfp', 'group-level/lme_prearousalsigma_placebo/input/')
end
if exist('group-level/lme_prearousalsigma_etc120', 'dir') == 0
    mkdir('group-level/lme_prearousalsigma_etc120/input')
    mkdir('group-level/lme_prearousalsigma_etc120/output')
    copyfile('group-level/paired-t-test_clustermass_perm/input/surface.srf', 'group-level/lme_prearousalsigma_etc120/input/')
    copyfile('group-level/paired-t-test_clustermass_perm/input/chanlocs.sfp', 'group-level/lme_prearousalsigma_etc120/input/')
end
% -------------------------------------------------------------------------
% Selection index
times = -90:0.5:30;
specfreqs = 1:1:30;
t_sel = [...
    -12.75, -4.25; ... % fast spindle
    -8.75, -4.24; ... % slow spindle
    ];
f_sel = [...
    14.5, 16.5; ... % fast spindle
    10.5, 14.5; ... % slow spindle
    ];
% -------------------------------------------------------------------------
% Create first level output files for PALM
for condition = {'placebo', 'etc120'}
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Select condition
    idx_cond_cs = strcmpi(Cond.cs, condition{:});
    idx_cond_aw = strcmpi(Cond.aw, condition{:});
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Fast spindle (depvar_1)
    idx_t = times > t_sel(1, 1) & times < t_sel(1, 2);
    idx_f = specfreqs > f_sel(1, 1) & specfreqs < f_sel(1, 2);
    depvar_1 = [...
        squeeze(mean(mean(CData.cs(:, idx_t, idx_f, idx_cond_cs), 2, 'omitnan'), 3, 'omitnan'))'; ...
        squeeze(mean(mean(CData.aw(:, idx_t, idx_f, idx_cond_aw), 2, 'omitnan'), 3, 'omitnan'))'; ...
        ];
    depvar_1 = log10(depvar_1);
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Slow spindle (depvar_2)
    idx_t = times > t_sel(2, 1) & times < t_sel(2, 2);
    idx_f = specfreqs > f_sel(2, 1) & specfreqs < f_sel(2, 2);
    depvar_2 = [...
        squeeze(mean(mean(CData.cs(:, idx_t, idx_f, idx_cond_cs), 2, 'omitnan'), 3, 'omitnan'))'; ...
        squeeze(mean(mean(CData.aw(:, idx_t, idx_f, idx_cond_aw), 2, 'omitnan'), 3, 'omitnan'))'; ...
        ];
    depvar_2 = log10(depvar_2);
    % -------------------------------------------------------------------------
    % Variance groups and exchangeability blocks
    vargroup = ones(sum(idx_cond_cs) + sum(idx_cond_aw), 1);
    exchblock = zeros(sum(idx_cond_cs) + sum(idx_cond_aw), 1);
    % -------------------------------------------------------------------------
    % Design matrix
    unique_subj = unique([Subj.cs; Subj.aw]);
    fix_fx = [-1 .* ones(sum(idx_cond_cs), 1); 1 .* ones(sum(idx_cond_aw), 1)];
    rand_fx = zeros(size(fix_fx, 1), length(unique_subj));
    for i = 1:length(unique_subj)
        rand_fx(strcmpi([Subj.cs(idx_cond_cs); Subj.aw(idx_cond_aw)], unique_subj{i}), i) = 1;
        exchblock(strcmpi([Subj.cs(idx_cond_cs); Subj.aw(idx_cond_aw)], unique_subj{i})) = i;
    end
    desmat = [ ...
        fix_fx, ...
        rand_fx, ...
        ];
    % -------------------------------------------------------------------------
    % Contrast
    tcontrast1 = zeros(1, size(desmat, 2));
    tcontrast1(1) = 1;
    % -------------------------------------------------------------------------
    % Save
    writematrix(depvar_1, sprintf('group-level/lme_prearousalsigma_%s/input/depvar_1.csv', condition{:}))
    writematrix(depvar_2, sprintf('group-level/lme_prearousalsigma_%s/input/depvar_2.csv', condition{:}))
    writematrix(desmat, sprintf('group-level/lme_prearousalsigma_%s/input/desmat.csv', condition{:}))
    writematrix(exchblock, sprintf('group-level/lme_prearousalsigma_%s/input/exchangeability_blocks.csv', condition{:}))
    writematrix(vargroup, sprintf('group-level/lme_prearousalsigma_%s/input/variance_groups.csv', condition{:}))
    writematrix(tcontrast1, sprintf('group-level/lme_prearousalsigma_%s/input/tcontrast1.csv', condition{:}))
    % Run
    cmd = ['palm -n 10000 ', ...
        sprintf('-i ''%s/group-level/lme_prearousalsigma_%s/input/depvar_1.csv'' ', pwd, condition{:}), ...
        sprintf('-i ''%s/group-level/lme_prearousalsigma_%s/input/depvar_2.csv'' ', pwd, condition{:}), ...
        sprintf('-s ''%s/group-level/lme_prearousalsigma_%s/input/surface.srf'' ', pwd, condition{:}), ...
        sprintf('-d ''%s/group-level/lme_prearousalsigma_%s/input/desmat.csv'' ', pwd, condition{:}), ...
        sprintf('-eb ''%s/group-level/lme_prearousalsigma_%s/input/exchangeability_blocks.csv'' ', pwd, condition{:}), ...
        sprintf('-vg ''%s/group-level/lme_prearousalsigma_%s/input/variance_groups.csv'' ', pwd, condition{:}), ...
        sprintf('-t ''%s/group-level/lme_prearousalsigma_%s/input/tcontrast1.csv'' ', pwd, condition{:}), ...
        sprintf('-o ''%s/group-level/lme_prearousalsigma_%s/output/palm_tstat1'' ', pwd, condition{:}), ...
        '-ee ', ...
        '-Cstat mass ', ...
        sprintf('-C %.4f ', norminv(1-0.001/2)), ...
        '-within ', ...
        '-twotail ', ...
        '-seed 251312519 ', ...
        '-cmcx ', ...
        '-savedof ', ...
        '-savemetrics ', ...
        '-saveglm ', ...
        '-verbosefilenames ', ...
        '-quiet'];
    eval(cmd)
end
% -------------------------------------------------------------------------
% Chanlocs and surface
chanlocs = template_to_chanlocs(which('GSN-HydroCel-257.sfp'));
[~, eegsystem] = ishdeeg({chanlocs.labels});
incl = hdeeg_scalpchannels(eegsystem);
chanlocs = chanlocs(ismember({chanlocs.labels}, incl));
surface = chanlocs2surface(...
    fullfile('group-level/lme_prearousalsigma_placebo/input', 'surface.srf'), ...
    fullfile('group-level/lme_prearousalsigma_placebo/input', 'chanlocs.sfp'), false, false);
%% -------------------------------------------------------------------------
% Load GLM results
PALM_NREM = struct();
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cfg = struct();
cfg.OutputDir = 'group-level/lme_prearousalsigma_placebo';
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
PALM_NREM.plc(1).label = 'fast-spindle';
[PALM_NREM.plc(1).tstat, PALM_NREM.plc(1).pval] = get_pvalue(cfg, 'T', 'p_unc', 'dpv', 1, 1);
[tstat, pval] = get_pvalue(cfg, 'T', 'p_fwe', 'clusterm', 1, 1);
PALM_NREM.plc(1).cluster = get_cluster(tstat, pval, 0.05, surface, chanlocs);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
PALM_NREM.plc(2).label = 'slow-spindle';
[PALM_NREM.plc(2).tstat, PALM_NREM.plc(2).pval] = get_pvalue(cfg, 'T', 'p_unc', 'dpv', 2, 1);
[tstat, pval] = get_pvalue(cfg, 'T', 'p_fwe', 'clusterm', 2, 1);
PALM_NREM.plc(2).cluster = get_cluster(tstat, pval, 0.05, surface, chanlocs);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cfg = struct();
cfg.OutputDir = 'group-level/lme_prearousalsigma_etc120';
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
PALM_NREM.etc(1).label = 'fast-spindle';
[PALM_NREM.etc(1).tstat, PALM_NREM.etc(1).pval] = get_pvalue(cfg, 'T', 'p_unc', 'dpv', 1, 1);
[tstat, pval] = get_pvalue(cfg, 'T', 'p_fwe', 'clusterm', 1, 1);
PALM_NREM.etc(1).cluster = get_cluster(tstat, pval, 0.05, surface, chanlocs);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
PALM_NREM.etc(2).label = 'slow-spindle';
[PALM_NREM.etc(2).tstat, PALM_NREM.etc(2).pval] = get_pvalue(cfg, 'T', 'p_unc', 'dpv', 2, 1);
[tstat, pval] = get_pvalue(cfg, 'T', 'p_fwe', 'clusterm', 2, 1);
PALM_NREM.etc(2).cluster = get_cluster(tstat, pval, 0.05, surface, chanlocs);

%%
% -------------------------------------------------------------------------
% REM SLEEP FINDINGS
% -------------------------------------------------------------------------
% Create dirs if needed
if exist('group-level/lme_prearousal_rem', 'dir') == 0
    mkdir('group-level/lme_prearousal_rem/input')
    mkdir('group-level/lme_prearousal_rem/output')
    copyfile('group-level/paired-t-test_clustermass_perm/input/surface.srf', 'group-level/lme_prearousal_rem/input/')
    copyfile('group-level/paired-t-test_clustermass_perm/input/chanlocs.sfp', 'group-level/lme_prearousal_rem/input/')
end
% -------------------------------------------------------------------------
% Selection index
times = -90:0.5:30;
specfreqs = 1:1:30;
t_sel = [...
    -18.75, -11.25; ... % beta
    -29.25, -9.75; ... % sigma
    -15.25, -8.5; ... % alpha
    ];
f_sel = [...
    15.5, 24.5; ... % beta
    11.5, 15.5; ... % sigma
    6.5, 11.5; ... % alpha
    ];
% -------------------------------------------------------------------------
% Create first level output files for PALM
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% beta (depvar_1)
idx_t = times > t_sel(1, 1) & times < t_sel(1, 2);
idx_f = specfreqs > f_sel(1, 1) & specfreqs < f_sel(1, 2);
depvar_1 = [...
    squeeze(mean(mean(CData.cs(:, idx_t, idx_f, :), 2, 'omitnan'), 3, 'omitnan'))'; ...
    squeeze(mean(mean(CData.aw(:, idx_t, idx_f, :), 2, 'omitnan'), 3, 'omitnan'))'; ...
    ];
depvar_1 = log10(depvar_1);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% sigma (depvar_2)
idx_t = times > t_sel(2, 1) & times < t_sel(2, 2);
idx_f = specfreqs > f_sel(2, 1) & specfreqs < f_sel(2, 2);
depvar_2 = [...
    squeeze(mean(mean(CData.cs(:, idx_t, idx_f, :), 2, 'omitnan'), 3, 'omitnan'))'; ...
    squeeze(mean(mean(CData.aw(:, idx_t, idx_f, :), 2, 'omitnan'), 3, 'omitnan'))'; ...
    ];
depvar_2 = log10(depvar_2);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% alpha (depvar_3)
idx_t = times > t_sel(3, 1) & times < t_sel(3, 2);
idx_f = specfreqs > f_sel(3, 1) & specfreqs < f_sel(3, 2);
depvar_3 = [...
    squeeze(mean(mean(CData.cs(:, idx_t, idx_f, :), 2, 'omitnan'), 3, 'omitnan'))'; ...
    squeeze(mean(mean(CData.aw(:, idx_t, idx_f, :), 2, 'omitnan'), 3, 'omitnan'))'; ...
    ];
depvar_3 = log10(depvar_3);
% -------------------------------------------------------------------------
% Variance groups and exchangeability blocks
vargroup = ones(size(depvar_1, 1), 1);
exchblock = zeros(size(depvar_1, 1), 1);
% -------------------------------------------------------------------------
% Design matrix
idx_cond_cs = strcmpi(Cond.cs, 'etc120');
idx_cond_aw = strcmpi(Cond.aw, 'etc120');
unique_subj = unique([Subj.cs; Subj.aw]);
fix_fx = [...
    round([idx_cond_cs; idx_cond_aw] - 0.5), ...% Main effect of treatment condition
    [-1 .* ones(length(idx_cond_cs), 1); 1 .* ones(length(idx_cond_aw), 1)], ... % Main effect of sleep stability
    ];
rand_fx = zeros(size(fix_fx, 1), length(unique_subj));
for i = 1:length(unique_subj)
    rand_fx(strcmpi([Subj.cs; Subj.aw], unique_subj{i}), i) = 1;
    exchblock(strcmpi([Subj.cs; Subj.aw], unique_subj{i})) = i;
end
desmat = [ ...
    fix_fx, ...
    rand_fx, ...
    ];
% -------------------------------------------------------------------------
% Contrast
tcontrast1 = zeros(1, size(desmat, 2));
tcontrast1(2) = 1; % second column in design matrix encodes sleep instability outcome
% -------------------------------------------------------------------------
% Save
writematrix(depvar_1, 'group-level/lme_prearousal_rem/input/depvar_1.csv')
writematrix(depvar_2, 'group-level/lme_prearousal_rem/input/depvar_2.csv')
writematrix(depvar_3, 'group-level/lme_prearousal_rem/input/depvar_3.csv')
writematrix(desmat, 'group-level/lme_prearousal_rem/input/desmat.csv')
writematrix(exchblock, 'group-level/lme_prearousal_rem/input/exchangeability_blocks.csv')
writematrix(vargroup, 'group-level/lme_prearousal_rem/input/variance_groups.csv')
writematrix(tcontrast1, 'group-level/lme_prearousal_rem/input/tcontrast1.csv')
% Run
cmd = ['palm -n 10000 ', ...
    sprintf('-i ''%s/group-level/lme_prearousal_rem/input/depvar_1.csv'' ', pwd), ...
    sprintf('-i ''%s/group-level/lme_prearousal_rem/input/depvar_2.csv'' ', pwd), ...
    sprintf('-i ''%s/group-level/lme_prearousal_rem/input/depvar_3.csv'' ', pwd), ...
    sprintf('-s ''%s/group-level/lme_prearousal_rem/input/surface.srf'' ', pwd), ...
    sprintf('-d ''%s/group-level/lme_prearousal_rem/input/desmat.csv'' ', pwd), ...
    sprintf('-eb ''%s/group-level/lme_prearousal_rem/input/exchangeability_blocks.csv'' ', pwd), ...
    sprintf('-vg ''%s/group-level/lme_prearousal_rem/input/variance_groups.csv'' ', pwd), ...
    sprintf('-t ''%s/group-level/lme_prearousal_rem/input/tcontrast1.csv'' ', pwd), ...
    sprintf('-o ''%s/group-level/lme_prearousal_rem/output/palm_tstat1'' ', pwd), ...
    '-ee ', ...
    '-Cstat mass ', ...
    sprintf('-C %.4f ', norminv(1-0.001/2)), ...
    '-within ', ...
    '-twotail ', ...
    '-seed 251312519 ', ...
    '-cmcx ', ...
    '-savedof ', ...
    '-savemetrics ', ...
    '-saveglm ', ...
    '-verbosefilenames ', ...
    '-quiet'];
eval(cmd)
% -------------------------------------------------------------------------
% Chanlocs and surface
chanlocs = template_to_chanlocs(which('GSN-HydroCel-257.sfp'));
[~, eegsystem] = ishdeeg({chanlocs.labels});
incl = hdeeg_scalpchannels(eegsystem);
chanlocs = chanlocs(ismember({chanlocs.labels}, incl));
surface = chanlocs2surface(...
    fullfile('group-level/lme_prearousal_rem/input', 'surface.srf'), ...
    fullfile('group-level/lme_prearousal_rem/input', 'chanlocs.sfp'), false, false);
% -------------------------------------------------------------------------
% Load GLM results
PALM_REM = struct();
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cfg = struct();
cfg.OutputDir = 'group-level/lme_prearousal_rem';
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
PALM_REM.main(1).label = 'beta';
[PALM_REM.main(1).tstat, PALM_REM.main(1).pval] = get_pvalue(cfg, 'T', 'p_unc', 'dpv', 1, 1);
[tstat, pval] = get_pvalue(cfg, 'T', 'p_fwe', 'clusterm', 1, 1);
PALM_REM.main(1).cluster = get_cluster(tstat, pval, 0.05, surface, chanlocs);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
PALM_REM.main(2).label = 'sigma';
[PALM_REM.main(2).tstat, PALM_REM.main(2).pval] = get_pvalue(cfg, 'T', 'p_unc', 'dpv', 2, 1);
[tstat, pval] = get_pvalue(cfg, 'T', 'p_fwe', 'clusterm', 2, 1);
PALM_REM.main(2).cluster = get_cluster(tstat, pval, 0.05, surface, chanlocs);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
PALM_REM.main(3).label = 'alpha';
[PALM_REM.main(3).tstat, PALM_REM.main(3).pval] = get_pvalue(cfg, 'T', 'p_unc', 'dpv', 3, 1);
[tstat, pval] = get_pvalue(cfg, 'T', 'p_fwe', 'clusterm', 3, 1);
PALM_REM.main(3).cluster = get_cluster(tstat, pval, 0.05, surface, chanlocs);

%% TOPOPLOTS FOR PALM RESULTS
% -------------------------------------------------------------------------
% NREM2 SLEEP FINDINGS
% -------------------------------------------------------------------------
% Selection index
times = -90:0.5:30;
specfreqs = 1:1:30;
t_sel = [...
    -12.75, -4.25; ... % fast spindle
    -8.75, -4.24; ... % slow spindle
    ];
f_sel = [...
    14.5, 16.5; ... % fast spindle
    10.5, 14.5; ... % slow spindle
    ];
% -------------------------------------------------------------------------
% PLACEBO
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Fig = figure();
Fig.Units = 'centimeters';
Fig.Position(3:4) = [10, 10];
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
CMapRoma = load('colormap_roma.mat');
CMapRoma = CMapRoma.roma;
CMapBatlow = load('colormap_batlow.mat');
CMapBatlow = CMapBatlow.batlow;
clear Ax
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Select condition
condition = 'placebo';
idx_cond_cs = strcmpi(Cond.cs, condition);
idx_cond_aw = strcmpi(Cond.aw, condition);
switch condition
    case 'placebo'
        fld = 'plc';
    case 'etc120'
        fld = 'etc';
end
% -------------------------------------------------------------------------
% Fast spindle
idx_t = times > t_sel(1, 1) & times < t_sel(1, 2);
idx_f = specfreqs > f_sel(1, 1) & specfreqs < f_sel(1, 2);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Slow spindle - continued sleep
Ax(1) = axes(Fig);
Ax(1).Position = [0.1 1.04-1/3 1/4 1/4];
Ax(1).YLabel.String = 'CONT. SLEEP';
YData = squeeze(median(sum(mean(CData.cs(:, idx_t, idx_f, idx_cond_cs), 2, 'omitnan'), 3, 'omitnan'), 4, 'omitnan'));
topoplot(log10(YData), chanlocs, ...
    'contourvals', (YData > 0), ...
    'numcontour', length(unique((YData > 0)))-1, ...
    'conv', 'on', ...
    'whitebk', 'on');
Ax(1).Title.String = 'FAST SIGMA';
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Slow spindle - state shift
Ax(2) = axes(Fig);
Ax(2).Position = [0.1 1-2/3 1/4 1/4];
Ax(2).YLabel.String = 'STATE SHIFT';
YData = squeeze(median(sum(mean(CData.aw(:, idx_t, idx_f, idx_cond_aw), 2, 'omitnan'), 3, 'omitnan'), 4, 'omitnan'));
topoplot(log10(YData), chanlocs, ...
    'contourvals', (YData > 0), ...
    'numcontour', length(unique((YData > 0)))-1, ...
    'conv', 'on', ...
    'whitebk', 'on');
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Slow spindle - t-stat
Ax(3) = axes(Fig);
Ax(3).Position = [0.1 1-3/3 1/4 1/4];
Ax(3).YLabel.String = 'T-STAT';
YData = PALM_NREM.(fld)(1).tstat;
if ~isempty(PALM_NREM.(fld)(1).cluster)
    clustlabels = arrayfun(@(c){c.chanlocs.labels}, PALM_NREM.(fld)(1).cluster, 'UniformOutput', false);
    clustlabels = unique([clustlabels{:}]);
    EMarkers = find(ismember({chanlocs.labels}, clustlabels));
else
    EMarkers = [];
end
topoplot(YData, chanlocs, ...
    'emarker2', {EMarkers, '.', 'w', 8, 1}, ... 
    'contourvals', (YData > -3.89) + (YData > -3.29) + (YData > -2.576) + (YData > -1.96) + (YData > 1.96) + (YData > 2.576) + (YData > 3.29) + (YData > 3.89), ...
    'numcontour', length(unique((YData > -3.89) + (YData > -3.29) + (YData > -2.576) + (YData > -1.96) + (YData > 1.96) + (YData > 2.576) + (YData > 3.29) + (YData > 3.89)))-1, ...
    'conv', 'on', ...
    'whitebk', 'on');
% -------------------------------------------------------------------------
% Slow spindle
idx_t = times > t_sel(2, 1) & times < t_sel(2, 2);
idx_f = specfreqs > f_sel(2, 1) & specfreqs < f_sel(2, 2);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Slow spindle - continued sleep
Ax(4) = axes(Fig);
Ax(4).Position = [0.4 1.04-1/3 1/4 1/4];
YData = squeeze(median(sum(mean(CData.cs(:, idx_t, idx_f, idx_cond_cs), 2, 'omitnan'), 3, 'omitnan'), 4, 'omitnan'));
topoplot(log10(YData), chanlocs, ...
    'contourvals', (YData > 0), ...
    'numcontour', length(unique((YData > 0)))-1, ...
    'conv', 'on', ...
    'whitebk', 'on');
tmp = YData;
Ax(4).Title.String = 'SLOW SIGMA';
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Slow spindle - state shift
Ax(5) = axes(Fig);
Ax(5).Position = [0.4 1-2/3 1/4 1/4];
YData = squeeze(median(sum(mean(CData.aw(:, idx_t, idx_f, idx_cond_aw), 2, 'omitnan'), 3, 'omitnan'), 4, 'omitnan'));
topoplot(log10(YData), chanlocs, ...
    'contourvals', (YData > 0), ...
    'numcontour', length(unique((YData > 0)))-1, ...
    'conv', 'on', ...
    'whitebk', 'on');
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Slow spindle - t-stat
Ax(6) = axes(Fig);
Ax(6).Position = [0.4 1-3/3 1/4 1/4];
YData = PALM_NREM.(fld)(2).tstat;
if ~isempty(PALM_NREM.(fld)(2).cluster)
    clustlabels = arrayfun(@(c){c.chanlocs.labels}, PALM_NREM.(fld)(2).cluster, 'UniformOutput', false);
    clustlabels = unique([clustlabels{:}]);
    EMarkers = find(ismember({chanlocs.labels}, clustlabels));
else
    EMarkers = [];
end
topoplot(YData, chanlocs, ...
    'emarker2', {EMarkers, '.', 'w', 8, 1}, ... 
    'contourvals', (YData > -3.89) + (YData > -3.29) + (YData > -2.576) + (YData > -1.96) + (YData > 1.96) + (YData > 2.576) + (YData > 3.29) + (YData > 3.89), ...
    'numcontour', length(unique((YData > -3.89) + (YData > -3.29) + (YData > -2.576) + (YData > -1.96) + (YData > 1.96) + (YData > 2.576) + (YData > 3.29) + (YData > 3.89)))-1, ...
    'conv', 'on', ...
    'whitebk', 'on');
% -------------------------------------------------------------------------
Ax(1).CLim = log10([0.0008 0.002]);
Ax(2).CLim = log10([0.0008 0.002]);
Ax(3).CLim = [-3.30 3.30];
Ax(4).CLim = log10([0.002 0.007]);
Ax(5).CLim = log10([0.002 0.007]);
Ax(6).CLim = [-3.30 3.30];
Ax(1).Colormap = CMapBatlow;
Ax(2).Colormap = CMapBatlow;
Ax(3).Colormap = CMapRoma;
Ax(4).Colormap = CMapBatlow;
Ax(5).Colormap = CMapBatlow;
Ax(6).Colormap = CMapRoma;
Ax(1).XColor = 'w';
Ax(2).XColor = 'w';
Ax(3).XColor = 'w';
Ax(1).YColor = 'w';
Ax(2).YColor = 'w';
Ax(3).YColor = 'w';
Ax(1).YLabel.Color = 'k';
Ax(2).YLabel.Color = 'k';
Ax(3).YLabel.Color = 'k';
Ax(1).YLabel.Position(1) = -0.67;
Ax(2).YLabel.Position(1) = -0.67;
Ax(3).YLabel.Position(1) = -0.67;
Ax(1).Visible = 'on';
Ax(2).Visible = 'on';
Ax(3).Visible = 'on';
% -------------------------------------------------------------------------
clear CBar
CBar(1) = colorbar(Ax(1), 'west');
CBar(1).Position = [0.21, 0.6, 0.03, 0.1];
CBar(1).Ticks = log10([0.0008 0.002]);
CBar(1).TickLabels = {'0.0008', '0.002'};
CBar(2) = colorbar(Ax(4), 'west');
CBar(2).Position = [0.51, 0.6, 0.03, 0.1];
CBar(2).Ticks = log10([0.002 0.007]);
CBar(2).TickLabels = {'0.002', '0.007'};
CBar(3) = colorbar(Ax(3), 'west');
CBar(3).Position = [0.7, 0.025, 0.03, 0.2];
CBar(3).Ticks = sort([-1*norminv(1-[0.05, 0.01, 0.001]/2), norminv(1-[0.05, 0.01, 0.001]/2)]);
CBar(3).TickLabels = arrayfun(@(v) sprintf('%.2f', v), CBar(3).Ticks, 'UniformOutput', false);
CBar(3).TickLength = 0.15;
drawnow()
Ax(1).Title.FontSize = 8;
Ax(4).Title.FontSize = 8;
Ax(1).YLabel.FontSize = 8;
Ax(2).YLabel.FontSize = 8;
Ax(3).YLabel.FontSize = 8;
CBar(1).FontSize = 6;
CBar(2).FontSize = 6;
CBar(3).FontSize = 6;

exportgraphics(Fig, sprintf('figures/fig_2_topoplots_nrem2_%s.png', condition), 'Resolution', 1200)

%%
% -------------------------------------------------------------------------
% REM SLEEP FINDINGS
% -------------------------------------------------------------------------
% Selection index
times = -90:0.5:30;
specfreqs = 1:1:30;
t_sel = [...
    -18.75, -11.25; ... % beta
    -29.25, -9.75; ... % sigma
    -15.25, -8.5; ... % alpha
    ];
f_sel = [...
    15.5, 24.5; ... % beta
    11.5, 15.5; ... % sigma
    6.5, 11.5; ... % alpha
    ];
% -------------------------------------------------------------------------
Fig = figure();
Fig.Units = 'centimeters';
Fig.Position(3:4) = [15, 10];
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
CMapRoma = load('colormap_roma.mat');
CMapRoma = CMapRoma.roma;
CMapHawaii = load('colormap_hawaii.mat');
CMapHawaii = CMapHawaii.hawaii;
CMapBatlow = load('colormap_batlow.mat');
CMapBatlow = CMapBatlow.batlow;
clear Ax
% -------------------------------------------------------------------------
% Beta
idx_t = times > t_sel(1, 1) & times < t_sel(1, 2);
idx_f = specfreqs > f_sel(1, 1) & specfreqs < f_sel(1, 2);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Beta - continued sleep
Ax(1) = axes(Fig);
Ax(1).Position = [0.1 1.04-1/3 1/4 1/4];
Ax(1).YLabel.String = 'CONT. SLEEP';
YData = squeeze(median(sum(mean(CData.cs(:, idx_t, idx_f, :), 2, 'omitnan'), 3, 'omitnan'), 4, 'omitnan'));
topoplot(log10(YData), chanlocs, ...
    'contourvals', (YData > 0), ...
    'numcontour', length(unique((YData > 0)))-1, ...
    'conv', 'on', ...
    'whitebk', 'on');
Ax(1).Title.String = 'BETA';
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Beta - state shift
Ax(2) = axes(Fig);
Ax(2).Position = [0.1 1-2/3 1/4 1/4];
Ax(2).YLabel.String = 'STATE SHIFT';
YData = squeeze(median(sum(mean(CData.aw(:, idx_t, idx_f, :), 2, 'omitnan'), 3, 'omitnan'), 4, 'omitnan'));
topoplot(log10(YData), chanlocs, ...
    'contourvals', (YData > 0), ...
    'numcontour', length(unique((YData > 0)))-1, ...
    'conv', 'on', ...
    'whitebk', 'on');
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Beta - t-stat
Ax(3) = axes(Fig);
Ax(3).Position = [0.1 1-3/3 1/4 1/4];
Ax(3).YLabel.String = 'T-STAT';
YData = PALM_REM.main(1).tstat;
if ~isempty(PALM_REM.main(1).cluster)
    clustlabels = arrayfun(@(c){c.chanlocs.labels}, PALM_REM.main(1).cluster, 'UniformOutput', false);
    clustlabels = unique([clustlabels{:}]);
    EMarkers = find(ismember({chanlocs.labels}, clustlabels));
else
    EMarkers = [];
end
topoplot(YData, chanlocs, ...
    'emarker2', {EMarkers, '.', 'w', 8, 1}, ... 
    'contourvals', (YData > -3.89) + (YData > -3.29) + (YData > -2.576) + (YData > -1.96) + (YData > 1.96) + (YData > 2.576) + (YData > 3.29) + (YData > 3.89), ...
    'numcontour', length(unique((YData > -3.89) + (YData > -3.29) + (YData > -2.576) + (YData > -1.96) + (YData > 1.96) + (YData > 2.576) + (YData > 3.29) + (YData > 3.89)))-1, ...
    'conv', 'on', ...
    'whitebk', 'on');
% -------------------------------------------------------------------------
% Sigma
idx_t = times > t_sel(2, 1) & times < t_sel(2, 2);
idx_f = specfreqs > f_sel(2, 1) & specfreqs < f_sel(2, 2);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Sigma - continued sleep
Ax(4) = axes(Fig);
Ax(4).Position = [0.3 1.04-1/3 1/4 1/4];
YData = squeeze(median(sum(mean(CData.cs(:, idx_t, idx_f, :), 2, 'omitnan'), 3, 'omitnan'), 4, 'omitnan'));
topoplot(log10(YData), chanlocs, ...
    'contourvals', (YData > 0), ...
    'numcontour', length(unique((YData > 0)))-1, ...
    'conv', 'on', ...
    'whitebk', 'on');
tmp = YData;
Ax(4).Title.String = 'SIGMA';
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Sigma - state shift
Ax(5) = axes(Fig);
Ax(5).Position = [0.3 1-2/3 1/4 1/4];
YData = squeeze(median(sum(mean(CData.aw(:, idx_t, idx_f, :), 2, 'omitnan'), 3, 'omitnan'), 4, 'omitnan'));
topoplot(log10(YData), chanlocs, ...
    'contourvals', (YData > 0), ...
    'numcontour', length(unique((YData > 0)))-1, ...
    'conv', 'on', ...
    'whitebk', 'on');
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Sigma - t-stat
Ax(6) = axes(Fig);
Ax(6).Position = [0.3 1-3/3 1/4 1/4];
YData = PALM_REM.main(2).tstat;
if ~isempty(PALM_REM.main(2).cluster)
    clustlabels = arrayfun(@(c){c.chanlocs.labels}, PALM_REM.main(2).cluster, 'UniformOutput', false);
    clustlabels = unique([clustlabels{:}]);
    EMarkers = find(ismember({chanlocs.labels}, clustlabels));
else
    EMarkers = [];
end
topoplot(YData, chanlocs, ...
    'emarker2', {EMarkers, '.', 'w', 8, 1}, ... 
    'contourvals', (YData > -3.89) + (YData > -3.29) + (YData > -2.576) + (YData > -1.96) + (YData > 1.96) + (YData > 2.576) + (YData > 3.29) + (YData > 3.89), ...
    'numcontour', length(unique((YData > -3.89) + (YData > -3.29) + (YData > -2.576) + (YData > -1.96) + (YData > 1.96) + (YData > 2.576) + (YData > 3.29) + (YData > 3.89)))-1, ...
    'conv', 'on', ...
    'whitebk', 'on');
% -------------------------------------------------------------------------
% Alpha
idx_t = times > t_sel(3, 1) & times < t_sel(3, 2);
idx_f = specfreqs > f_sel(3, 1) & specfreqs < f_sel(3, 2);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Alpha - continued sleep
Ax(7) = axes(Fig);
Ax(7).Position = [0.5 1.04-1/3 1/4 1/4];
YData = squeeze(median(sum(mean(CData.cs(:, idx_t, idx_f, :), 2, 'omitnan'), 3, 'omitnan'), 4, 'omitnan'));
topoplot(log10(YData), chanlocs, ...
    'contourvals', (YData > 0), ...
    'numcontour', length(unique((YData > 0)))-1, ...
    'conv', 'on', ...
    'whitebk', 'on');
tmp = YData;
Ax(7).Title.String = 'ALPHA';
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Alpha - state shift
Ax(8) = axes(Fig);
Ax(8).Position = [0.5 1-2/3 1/4 1/4];
YData = squeeze(median(sum(mean(CData.aw(:, idx_t, idx_f, :), 2, 'omitnan'), 3, 'omitnan'), 4, 'omitnan'));
topoplot(log10(YData), chanlocs, ...
    'contourvals', (YData > 0), ...
    'numcontour', length(unique((YData > 0)))-1, ...
    'conv', 'on', ...
    'whitebk', 'on');
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Alpha - t-stat
Ax(9) = axes(Fig);
Ax(9).Position = [0.5 1-3/3 1/4 1/4];
YData = PALM_REM.main(3).tstat;
if ~isempty(PALM_REM.main(3).cluster)
    clustlabels = arrayfun(@(c){c.chanlocs.labels}, PALM_REM.main(3).cluster, 'UniformOutput', false);
    clustlabels = unique([clustlabels{:}]);
    EMarkers = find(ismember({chanlocs.labels}, clustlabels));
else
    EMarkers = [];
end
topoplot(YData, chanlocs, ...
    'emarker2', {EMarkers, '.', 'w', 8, 1}, ... 
    'contourvals', (YData > -3.89) + (YData > -3.29) + (YData > -2.576) + (YData > -1.96) + (YData > 1.96) + (YData > 2.576) + (YData > 3.29) + (YData > 3.89), ...
    'numcontour', length(unique((YData > -3.89) + (YData > -3.29) + (YData > -2.576) + (YData > -1.96) + (YData > 1.96) + (YData > 2.576) + (YData > 3.29) + (YData > 3.89)))-1, ...
    'conv', 'on', ...
    'whitebk', 'on');

% -------------------------------------------------------------------------
Ax(1).CLim = log10([0.001 0.0025]);
Ax(2).CLim = log10([0.001 0.0025]);
Ax(3).CLim = [-3.30 3.30];
Ax(4).CLim = log10([0.001 0.004]);
Ax(5).CLim = log10([0.001 0.004]);
Ax(6).CLim = [-3.30 3.30];
Ax(7).CLim = log10([0.002 0.01]);
Ax(8).CLim = log10([0.002 0.01]);
Ax(9).CLim = [-3.30 3.30];
Ax(1).Colormap = CMapBatlow;
Ax(2).Colormap = CMapBatlow;
Ax(3).Colormap = CMapRoma;
Ax(4).Colormap = CMapBatlow;
Ax(5).Colormap = CMapBatlow;
Ax(6).Colormap = CMapRoma;
Ax(7).Colormap = CMapBatlow;
Ax(8).Colormap = CMapBatlow;
Ax(9).Colormap = CMapRoma;
Ax(1).XColor = 'w';
Ax(2).XColor = 'w';
Ax(3).XColor = 'w';
Ax(1).YColor = 'w';
Ax(2).YColor = 'w';
Ax(3).YColor = 'w';
Ax(1).YLabel.Color = 'k';
Ax(2).YLabel.Color = 'k';
Ax(3).YLabel.Color = 'k';
Ax(1).YLabel.Position(1) = -0.67;
Ax(2).YLabel.Position(1) = -0.67;
Ax(3).YLabel.Position(1) = -0.67;
Ax(1).Visible = 'on';
Ax(2).Visible = 'on';
Ax(3).Visible = 'on';
% -------------------------------------------------------------------------
clear CBar
CBar(1) = colorbar(Ax(1), 'west');
CBar(1).Position = [0.215, 0.6, 0.02, 0.1];
CBar(1).Ticks = log10([0.001 0.0025]);
CBar(1).TickLabels = {'0.001', '0.0025'};
CBar(2) = colorbar(Ax(4), 'west');
CBar(2).Position = [0.415, 0.6, 0.02, 0.1];
CBar(2).Ticks = log10([0.001 0.004]);
CBar(2).TickLabels = {'0.001', '0.004'};
CBar(3) = colorbar(Ax(7), 'west');
CBar(3).Position = [0.615, 0.6, 0.02, 0.1];
CBar(3).Ticks = log10([0.002 0.01]);
CBar(3).TickLabels = {'0.002', '0.01'};
CBar(4) = colorbar(Ax(3), 'west');
CBar(4).Position = [0.73, 0.025, 0.02, 0.2];
CBar(4).Ticks = sort([-1*norminv(1-[0.05, 0.01, 0.001]/2), norminv(1-[0.05, 0.01, 0.001]/2)]);
CBar(4).TickLabels = arrayfun(@(v) sprintf('%.2f', v), CBar(4).Ticks, 'UniformOutput', false);
CBar(4).TickLength = 0.15;
drawnow();
Ax(1).Title.FontSize = 8;
Ax(4).Title.FontSize = 8;
Ax(1).YLabel.FontSize = 8;
Ax(2).YLabel.FontSize = 8;
Ax(3).YLabel.FontSize = 8;
CBar(1).FontSize = 6;
CBar(2).FontSize = 6;
CBar(3).FontSize = 6;

exportgraphics(Fig, 'figures/fig_2_topoplots_rem.png', 'Resolution', 1200)

%% PRE-AROUSAL SPECTRUM
% -------------------------------------------------------------------------
% Parameters
times = -90:0.5:30;
specfreqs = 1:1:30;
XLim = [-30, 10];
CLim = [-4, -1];
RatioCLim = [-0.25, 0.25];
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Change here for the different conditions
sleepstage = 'rem';
fld_cond = 'plc';
if strcmpi(sleepstage, 'nrem2')
    switch fld_cond
        case 'plc'
            idx_cond_cs = strcmpi(Cond.cs, 'placebo');
            idx_cond_aw = strcmpi(Cond.aw, 'placebo');
        case 'etc'
            idx_cond_cs = ~strcmpi(Cond.cs, 'placebo');
            idx_cond_aw = ~strcmpi(Cond.aw, 'placebo');
    end
else
    idx_cond_cs = strcmpi(Cond.cs, 'placebo') | ~strcmpi(Cond.cs, 'placebo');
    idx_cond_aw = strcmpi(Cond.aw, 'placebo') | ~strcmpi(Cond.aw, 'placebo');
end
% -------------------------------------------------------------------------
% Create figure
Fig = figure();
Fig.Units = 'centimeters';
Fig.Position = [5, 5, 9, 10];
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
CMapRoma = load('colormap_roma.mat');
CMapRoma = CMapRoma.roma;
% -------------------------------------------------------------------------
% Percentage of each timepoint with arousals
Ax = axes(Fig);
Ax.NextPlot = 'add';
Ax.TickDir = 'out';
Ax.Layer = 'top';
Ax.Box = 'on';
Ax.Position = [0.3 0.73 0.55 0.1];
Ax.XLim = XLim;
Ax.XTick = sort([XLim, 0]);
Ax.XTickLabel = {''};
Ax.FontSize = 8;
Ax.YLabel.String = 'p_{arousal}';
Ax.YLabel.FontSize = 8;
Ax.YLim = [0, 1];
Ax.Title.FontSize = 8;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
patch(...
    'XData', [times(1), times, times(end)], ...
    'YData', [0, mean(AData.cs(idx_cond_cs, :), 1, 'omitnan'), 0], ...
    'LineStyle', 'none', ...
    'FaceColor', standard_colors('blue'), ...
    'FaceAlpha', 0.25)
patch(...
    'XData', [times(1), times, times(end)], ...
    'YData', [0, mean(AData.aw(idx_cond_aw, :), 1, 'omitnan'), 0], ...
    'LineStyle', 'none', ...
    'FaceColor', standard_colors('red'), ...
    'FaceAlpha', 0.25)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
YTick = max([mean(AData.cs(:, times < 0), 1, 'omitnan'), mean(AData.aw(:, times < 0), 1, 'omitnan')]);
Ax.YTick = [ceil(YTick*100)/100, 1];

% -------------------------------------------------------------------------
% Pre-arousal spectrum
Ax = axes(Fig);
Ax.NextPlot = 'add';
Ax.TickDir = 'out';
Ax.Layer = 'top';
Ax.Box = 'on';
Ax.XLim = 10.^CLim;
Ax.XTick = 10.^(CLim(1):CLim(2));
Ax.XTickLabel = [num2cell(Ax.XTick(1:end-1)), {''}];
Ax.XTickLabelRotation = 90;
Ax.XScale = 'log';
Ax.YLim = [0, 30];
Ax.YTick = [1, 4, 8, 12, 15, 25, 30];
Ax.YGrid = 'on';
Ax.Position = [0.1 0.3 0.2 0.4];
Ax.FontSize = 8;
Ax.YLabel.String = 'Frequency (Hz)';
Ax.YLabel.FontSize = 8;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
idx_t = times > XLim(1) & times < 0;
YData = specfreqs;
XData = mean(squeeze(mean(mean(CData.cs(:, idx_t, :, idx_cond_cs), 2, 'omitnan'), 4, 'omitnan')), 1, 'omitnan');
plot(Ax, XData, YData, ...
    'Color', standard_colors('blue'), ...
    'LineWidth', 2);
XData = mean(squeeze(mean(mean(CData.aw(:, idx_t, :, idx_cond_aw), 2, 'omitnan'), 4, 'omitnan')), 1, 'omitnan');
plot(Ax, XData, YData, ...
    'Color', standard_colors('red'), ...
    'LineWidth', 2);

% -------------------------------------------------------------------------
% Main spectrogram
Ax = axes(Fig);
Ax.NextPlot = 'add';
Ax.TickDir = 'out';
Ax.Layer = 'top';
Ax.Box = 'on';
Ax.XLim = XLim;
Ax.XTick = sort([XLim, 0]);
Ax.YLim = [0, 30];
Ax.CLim = RatioCLim;
Ax.XTick = [-90, 0, 30];
Ax.YTick = [1, 4, 8, 12, 15, 25, 30];
Ax.XTickLabel = {''};
Ax.YTickLabel = {''};
Ax.YGrid = 'on';
Ax.Position = [0.3 0.3 0.55 0.4];
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Extract the pre-arousal spectra of the two state-shift conditions
AvCDataCS = squeeze(median(CData.cs(:, :, :, idx_cond_cs), 4, 'omitnan')); % average across trials
AvCDataAW = squeeze(median(CData.aw(:, :, :, idx_cond_aw), 4, 'omitnan')); % average across trials
AvCData = log10(AvCDataAW ./ AvCDataCS);
AvCData = squeeze(mean(AvCData, 1, 'omitnan')); % average across channels
AvCData = AvCData'; % transpose
imagesc('XData', times, 'YData', specfreqs, 'CData', AvCData)
plot([0, 0], [0, 31], '-w')
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
XData = -90:0.5:-0.5;
YData = 1:30;
n_tests = 60*30; % number of correlated tests (60 timepoints and 30 freqs)
alpha_sidak = (1-normcdf(3.1))*2; %1-(1-0.05).^(1/n_tests);
alpha_posthoc = (1-normcdf(1.96)); %1-(1-0.05).^(1/n_tests);
switch sleepstage
    case 'nrem2'
        MainSData = double(squeeze(sum(GLM.logpow.awake.pval < alpha_sidak, 1)) > 0)';
        contour(XData, YData, MainSData, 1, 'k');
        IntxSData = double(squeeze(sum(GLM.logpow.intx.pval < alpha_sidak, 1)) > 0)';
        contour(XData, YData, IntxSData, 1, 'w', 'LineWidth', 1.5);
        PostHocSData = double(squeeze(sum(PostHoc.(fld_cond).awake.pval < alpha_posthoc, 1)) > 0)';
        PostHocSData = PostHocSData.*IntxSData;
        contour(XData, YData, PostHocSData, 1, 'k');
    case 'rem'
        IntxSData = double(squeeze(sum(GLM.logpow.intx.pval < alpha_sidak, 1)) > 0)';
        contour(XData, YData, IntxSData, 1, 'w', 'LineWidth', 1.5);
        MainSData = double(squeeze(sum(GLM.rem.awake.pval < alpha_sidak, 1)) > 0)';
        contour(XData, YData, MainSData, 1, 'k');
end
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
switch sleepstage
    case 'nrem2'
        plot(Ax, [-12.75, -4.25, -4.25, -12.75, -12.75], [14.5, 14.5, 16.5, 16.5, 14.5], '-', 'Color', standard_colors('green'))
        plot(Ax, [-8.75, -4.25, -4.25, -8.75, -8.75], [10.5, 10.5, 14.5, 14.5, 10.5], '-', 'Color', standard_colors('green'))
    case 'rem'
        plot(Ax, [-18.5, -11.25, -11.25, -18.5, -18.5], [15.5, 15.5, 24.5, 24.5, 15.5], '-', 'Color', standard_colors('green'))
        plot(Ax, [-29, -9.75, -9.75, -29, -29], [11.5, 11.5, 15.5, 15.5, 11.5], '-', 'Color', standard_colors('green'))
        plot(Ax, [-15.25, -8.5, -8.5, -15.25, -15.25], [6.5, 6.5, 11.5, 11.5, 6.5], '-', 'Color', standard_colors('green'))
end
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Ax.Colormap = CMapRoma;

% -------------------------------------------------------------------------
% Total power over time
Ax = axes(Fig);
Ax.NextPlot = 'add';
Ax.TickDir = 'out';
Ax.Layer = 'top';
Ax.Box = 'on';
Ax.XLim = XLim;
Ax.XTick = sort([XLim, 0]);
Ax.YLim = 10.^CLim;
Ax.YTick = 10.^(CLim(1):CLim(2));
Ax.YTickLabel = [num2cell(Ax.YTick(1:end-1)), {''}];
Ax.YScale = 'log';
Ax.XGrid = 'on';
Ax.Position = [0.3 0.15 0.55 0.15];
Ax.YAxisLocation = 'right';
Ax.FontSize = 8;
Ax.XLabel.String = 'Time (s)';
Ax.XLabel.FontSize = 8;
Ax.XLabel.Position(2) = 0.000045;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
XData = times;
YData = squeeze(mean(mean(mean(CData.cs(:, :, :, idx_cond_cs), 1, 'omitnan'), 3, 'omitnan'), 4, 'omitnan'));
p_leg(1) = plot(Ax, XData, YData, ...
    'Color', standard_colors('blue'), ...
    'LineWidth', 2);
YData = squeeze(mean(mean(mean(CData.aw(:, :, :, idx_cond_aw), 1, 'omitnan'), 3, 'omitnan'), 4, 'omitnan'));
p_leg(2) = plot(Ax, XData, YData, ...
    'Color', standard_colors('red'), ...
    'LineWidth', 2);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
leg_labels = {...
    sprintf('continued sleep (N = %i)', size(CData.cs(:, :, :, idx_cond_cs), 4)), ...
    sprintf('state shift (N = %i)', size(CData.aw(:, :, :, idx_cond_aw), 4)), ...
    };
leg = legend(Ax, p_leg, leg_labels, ...
    'Box', 'off', ...
    'FontSize', 8);
leg.Position(1) = Ax.Position(1) + Ax.Position(3) - leg.Position(3);
leg.Position(2) = 0.83;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
CBar = colorbar();
CBar.Location = 'south';
CBar.Position = [0.3, 0.05 0.55 0.03];
CBar.Limits = [0, 1];
CBar.Ticks = [0, 1];
CBar.TickLabels = RatioCLim;
CBar.Parent.Colormap = CMapRoma;

switch sleepstage
    case 'nrem2'
        figoutname = sprintf('figures/fig_2_%s_%s.png', sleepstage, fld_cond);
    case 'rem'
        figoutname = sprintf('figures/fig_2_%s.png', sleepstage);
end
exportgraphics(Fig, figoutname, 'Resolution', 1200)

