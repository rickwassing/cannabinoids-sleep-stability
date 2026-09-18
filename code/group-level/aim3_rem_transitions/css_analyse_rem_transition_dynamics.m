function css_analyse_rem_transition_dynamics()
% -------------------------------------------------------------------------
% Load all REM onsets and plot the ~100 seconds of sigma ISF
SigmaFiles = dir('derivatives/EEG-segmented/sub-*/ses-*/sub-*-sigmaprerembout_pow.set');
% -------------------------------------------------------------------------
% Settings
clear chans
chans.sel = 'pz';
chans.locs = template_to_chanlocs(which('GSN-HydroCel-257.sfp'));
chans.locs = chans.locs(ismember({chans.locs.labels}, hdeeg_scalpchannels('egi257')));
chans.locs = channel_clusters(chans.locs, 'mff');
chans.idx_p = [77, 78, 79, 85, 86, 87, 92, 93, 94, 99, 104, 110, 111, 112, 113, 120, 121, 122];
chans.idx_f = [21, 28, 13, 29, 22, 14, 5, 30, 23, 15, 6, 172, 24, 16, 7, 166, 17, 8, 160];
if strcmpi(chans.sel, 'fz')
    chans.idx = chans.idx_f;
else
    chans.idx = chans.idx_p;
end
% -------------------------------------------------------------------------
% Set filter config
SIGMA = LoadDataset(fullfile(SigmaFiles(1).folder, SigmaFiles(1).name), 'all');
% Emperical values for the ISF filter cutoffs
clear filtcfg
filtcfg.ford_mult = 2;
filtcfg.cutoff = [0.0138 0.0259];
filtcfg.wintype = 'kaiser';
filtcfg.transbw = 1/50;
filtcfg.rippledev = 0.05;
filtcfg.warg = 6;
filtcfg.order = pop_firwsord(filtcfg.wintype, SIGMA(1).srate, filtcfg.transbw, filtcfg.rippledev);
% -------------------------------------------------------------------------
% Settings for detecting troughs and peaks in the sigma ISF
evcfg = struct();
evcfg.srate = SIGMA(1).srate;
evcfg.pkdist = 30;
evcfg.pkheight = 0.75.*pi;
evcfg.mindelay = 20; % seconds. Minimum time prior to REM onset to identify the ISF trough
% -------------------------------------------------------------------------
% Load, filter, and extract sigma ISF peaks/troughs for every subject/bout
H = extract_prerem_isf_traces(SigmaFiles, chans, filtcfg, evcfg);

%% For each trace extract the peak amplitudes and duration
[avsigma, desmat] = build_rem_transition_features(H);

%%

clc

fprintf('Mean (SD) ISF half-duration in NREM: %.2f (%.2f) s.\n', mean(desmat.duration(strcmpi(desmat.stage, 'nrem'))), std(desmat.duration(strcmpi(desmat.stage, 'nrem'))))
fprintf('Mean (SD) ISF half-duration in TR: %.2f (%.2f) s.\n', mean(desmat.duration(strcmpi(desmat.stage, 'rem'))), std(desmat.duration(strcmpi(desmat.stage, 'rem'))))
fprintf('Mean (SD) ISF amplitude in Placebo: %.2f (%.2f) dB.\n', mean(desmat.amplitude(strcmpi(desmat.cond, 'placebo'))), std(desmat.amplitude(strcmpi(desmat.cond, 'placebo'))))
fprintf('Mean (SD) ISF amplitude in THC/CBD: %.2f (%.2f) dB.\n', mean(desmat.amplitude(strcmpi(desmat.cond, 'etc120'))), std(desmat.amplitude(strcmpi(desmat.cond, 'etc120'))))

%%
% Fit a linear model to check for differences in Sigma power prior to REM
% transitions

for it = 1:30
    [perms] = transrem_permutation_models_optimized(avsigma, H(1).times)
    save(sprintf('analysis_3_neurob_%s.mat', datestr(now, 'yyyymmddTHHMM')), 'perms', '-v7.3') %#ok<TNOW1,DATST>
end

%% Load all permutations and append
perm_files = dir('analysis_3*.mat');
clear perms
tr = now(); %#ok<TNOW1>
for i = 1:length(perm_files)
    tmp = load(fullfile(perm_files(i).folder, perm_files(i).name), 'perms');
    tr = remainingTime(tr, length(perm_files));
    if i == 1
        perms = tmp.perms;
        continue
    end
    perms = [perms; tmp.perms(2:end)]; %#ok<AGROW>
end

%% Crop
perms = perms(1:6000);

clear i tr tmp

%% Calculate cluster p-values (length of consecutive p-values < 0.05)

Clust = getpermpvalue(perms, 'y', '');

%%

clc
mdl = fitlme(desmat, 'amp_z ~ 1 + cond*stage + (1|sub)') %#ok<*NOPTS>

%

clc
mdl = fitlme(desmat, 'dur_z ~ 1 + cond*stage + (1|sub)')

%
clc
mdl = fitlme(desmat, 'remlat_z ~ 1 + amp_z*stage + cond + (1|sub)')

%%
clc
mdl = fitlme(desmat, 'remlat_z ~ 1 + dur_z*stage + cond + (1|sub)')
mdl_nrem = fitlme(desmat(strcmpi(desmat.stage, 'nrem'), :), 'remlat_z ~ 1 + dur_z + cond + (1|sub)')
mdl_rem = fitlme(desmat(strcmpi(desmat.stage, 'rem'), :), 'remlat_z ~ 1 + dur_z + cond + (1|sub)')

% mdl = fitlme(desmat(strcmpi(desmat.stage, 'nrem'), :), 'remlat_z ~ 1 + dur_z + (1|sub)')
% mdl = fitlme(desmat(strcmpi(desmat.stage, 'rem'), :), 'remlat_z ~ 1 + dur_z + (1|sub)')

%%

SIG = LoadDataset('./derivatives/EEG-processed/sub-r011/ses-placebo/sub-r011_ses-placebo_task-psg_desc-sigma_pow.set', 'header');
SIG.event(contains({SIG.event.type}, 'arousal') | contains({SIG.event.type}, 'alpha')) = [];
HYP = css_eeglab2hypnogram(SIG);
HYP.episode(HYP.episode == -1) = 1;
HYP.episode(HYP.episode == -2) = 1;
bouts = HYP.times(find(diff(double(HYP.episode == 2)) == 1));
bouts = [bouts-300/(60*60*24), bouts+60/(60*60*24)];

%%
% Assemble and export Figure 3
plot_fig3_full(SIG, HYP, bouts, chans, H, desmat, chans.sel);

end

