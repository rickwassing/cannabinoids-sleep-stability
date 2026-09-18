function css_analyse_arousal_isf_phase()
% -------------------------------------------------------------------------
% Use the the 130-second pre-arousal bouts of continuous N2 sleep to
% determine the phase angle, amplitude (hilbert) between state-shift
% arousals and continued sleep arousals; and whether CBN modulates this
% - State-shift arousals occur after the sigma ISF peak, and
%   continued-sleep arousals occur prior to the peak.
% - State-shift arousals occur at higher ISF amplitudes (hilbert)
% - CBN increases the number of arousals to after the peak, or at higher
%   ISF amplitudes
% -------------------------------------------------------------------------
% INIT
% -------------------------------------------------------------------------
clc
clear ALLCOM ALLEEG CURRENTSET CURRENTSTUDY EEG eeglabUpdater LASTCOM PLUGINLIST STUDY
close all
%% -------------------------------------------------------------------------
% Channel constants
% -------------------------------------------------------------------------
% Define the set of parietal channels to focus on
clear chans
chans.locs = template_to_chanlocs(which('GSN-HydroCel-257.sfp'));
chans.locs = chans.locs(ismember({chans.locs.labels}, hdeeg_scalpchannels('egi257')));
chans.locs = channel_clusters(chans.locs, 'mff');
chans.idx_p = [77, 78, 79, 85, 86, 87, 92, 93, 94, 99, 104, 110, 111, 112, 113, 120, 121, 122];
chans.idx_f = [21, 28, 13, 29, 22, 14, 5, 30, 23, 15, 6, 172, 24, 16, 7, 166, 17, 8, 160];

%% -------------------------------------------------------------------------
% Load segmented pre-arousal bouts
Files = dir('derivatives/EEG-segmented/sub-*/ses-*/sub-*_desc-sigmanremarobout*.set');
clear SIGMA SIGMA_P SIGMA_F
for i = 1:length(Files)
    SIGMA = LoadDataset(fullfile(Files(i).folder, Files(i).name), 'all');
    SIGMA_P(i) = pop_select(SIGMA, 'channel', chans.idx_p);
    SIGMA_F(i) = pop_select(SIGMA, 'channel', chans.idx_f);
end
clear i SIGMA Files

%% -------------------------------------------------------------------------
% Analysis constants
% -------------------------------------------------------------------------
% Crop the signal by 31.5 seconds from the end-boundary (30 seconds were
% added to include the arousal, and the wavelett at 11 Hz had a duration of
% 1.45 seconds (16 cycles * 1/11 seconds per cycle).
clear const
const.crop_aro = [0, -31.5].*SIGMA_P(1).srate;
const.append_type = 'zeros';
const.selaro = readtable('./inspect/selectedaros.csv', 'Delimiter', ',');
% -------------------------------------------------------------------------
% Emperical values for the ISF filter cutoffs
clear filtcfg
filtcfg.ford_mult = 2;
filtcfg.cutoff = [0.0138 0.0259];
filtcfg.wintype = 'kaiser';
filtcfg.transbw = 1/50;
filtcfg.rippledev = 0.05;
filtcfg.warg = 6;
filtcfg.order = pop_firwsord(filtcfg.wintype, SIGMA_P(1).srate, filtcfg.transbw, filtcfg.rippledev);

close all
[Fig, fname] = plot_filt_params(SIGMA_P, filtcfg, const);
exportgraphics(Fig, fname, 'Resolution', 300)

%% -------------------------------------------------------------------------
% Load EEG for drawing figure
Files = dir('./derivatives/EEG-preproc/sub-r019/ses-etc120/sub-r019_ses-etc120_task-psg_desc-plot_eeg.set');
EEG = LoadDataset(fullfile(Files(1).folder, Files(1).name), 'all');
SIG = ft_wavelettransform(EEG, waveletsettings(EEG));
SIG.data(isnan(SIG.data)) = 0;
FSIG = pop_firws(SIG, ...
    'fcutoff', filtcfg.cutoff, ...
    'ftype', 'bandpass', ...
    'wtype', filtcfg.wintype, ...
    'warg', filtcfg.warg, ...
    'forder', filtcfg.order, ...
    'plotfresp', false, ...
    'minphase', 0);
clear Files

%% -------------------------------------------------------------------------
% Init output table
T_pz = getarousalphaseangle(SIGMA_P, const, filtcfg); %#ok<NASGU>
T_fz = getarousalphaseangle(SIGMA_F, const, filtcfg); %#ok<NASGU>

%% or
load('./analysis_2b_20251118T2104.mat', 'T_fz', 'T_pz')

%%
% For each timepoint in rawdata test for differences using permutation tests
for it = 1:12
    [perms_pz] = prearousal_permutation_models_optimized(T_pz);
    [perms_fz] = prearousal_permutation_models_optimized(T_fz);

    save(sprintf('analysis_2b_neurob_%s.mat', datestr(now, 'yyyymmddTHHMM')), 'perms_pz', 'perms_fz', '-v7.3') %#ok<TNOW1,DATST>
end
disp('done')

%% Load all permutations and append
perm_files = dir('analysis_2b_*.mat');
clear perms_pz perms_fz
tr = now(); %#ok<TNOW1>
for i = 1:length(perm_files)
    tmp = load(fullfile(perm_files(i).folder, perm_files(i).name), 'perms_fz', 'perms_pz');
    tr = remainingTime(tr, length(perm_files));
    if i == 1
        perms_pz = tmp.perms_pz;
        perms_fz = tmp.perms_fz;
        continue
    end
    perms_pz = [perms_pz; tmp.perms_pz(2:end)]; %#ok<AGROW>
    perms_fz = [perms_fz; tmp.perms_fz(2:end)]; %#ok<AGROW>
end

% Crop
perms_pz = perms_pz(1:3000);
perms_fz = perms_fz(1:3000);

clear i tr tmp

%% Calculate cluster p-values (length of consecutive p-values < 0.05)
Clust_pz = struct();
Clust_pz.y.cs = getpermpvalue(perms_pz, 'y', 'cs');
Clust_pz.y.aw = getpermpvalue(perms_pz, 'y', 'aw');
Clust_pz.pr.cs = getpermpvalue(perms_pz, 'pr', 'cs');
Clust_pz.pr.aw = getpermpvalue(perms_pz, 'pr', 'aw');

Clust_fz = struct();
Clust_fz.y.cs = getpermpvalue(perms_fz, 'y', 'cs');
Clust_fz.y.aw = getpermpvalue(perms_fz, 'y', 'aw');
Clust_fz.pr.cs = getpermpvalue(perms_fz, 'pr', 'cs');
Clust_fz.pr.aw = getpermpvalue(perms_fz, 'pr', 'aw');

%%
clc
ts = -61.5:0.1:-1.5;
TableS1a = table();
TableS1a.Cluster = ascolumn(1:numel(Clust_fz.pr.aw.length));
TableS1a.Size = ascolumn(Clust_fz.pr.aw.length);
TableS1a.Time = ascolumn(arrayfun(@(s, e) sprintf('%.1f to %.1f s', ts(s), ts(e)), Clust_fz.pr.aw.idx(1,:), Clust_fz.pr.aw.idx(2,:), 'UniformOutput', false));
TableS1a.P = ascolumn(arrayfun(@(p) sprintf('%.3f', p), Clust_fz.pr.aw.perm_p, 'UniformOutput', false));
disp(TableS1a);

TableS1b = table();
TableS1b.Cluster = ascolumn(1:numel(Clust_fz.pr.cs.length));
TableS1b.Size = ascolumn(Clust_fz.pr.cs.length);
TableS1b.Time = ascolumn(arrayfun(@(s, e) sprintf('%.1f to %.1f s', ts(s), ts(e)), Clust_fz.pr.cs.idx(1,:), Clust_fz.pr.cs.idx(2,:), 'UniformOutput', false));
TableS1b.P = ascolumn(arrayfun(@(p) sprintf('%.3f', p), Clust_fz.pr.cs.perm_p, 'UniformOutput', false));
disp(TableS1b);

%% -------------------------------------------------------------------------
% Plot figure 2
% -------------------------------------------------------------------------
% Indexes
clc
clear pcfg
rois = {'fz', 'pz'};
delay = 'phase_d0';

Circ = table();
Circ.Y = {'AW-PBO'; 'AW-THC'; 'CS-PBO'; 'CS-THC'};
Circ.Ray_fz = cell(4,1);
Circ.HR_fz = cell(4,1);
Circ.Dip_fz = cell(4,1);
Circ.Ray_pz = cell(4,1);
Circ.HR_pz = cell(4,1);
Circ.Dip_pz = cell(4,1);

for ri = 1:length(rois)
    roi = rois{ri};

    fprintf('---------------------\n');
    fprintf('%s\n', roi);

    switch roi
        case 'fz'
            T_this = T_fz;
            Chans_this = chans.idx_f;
            Perms_this = perms_fz;
            Clust_this = Clust_fz;
        case 'pz'
            T_this = T_pz;
            Chans_this = chans.idx_p;
            Perms_this = perms_pz;
            Clust_this = Clust_pz;
    end

    T_this.(delay) = correct_phase_by_empirical_cdf(T_this.(delay));

    pcfg.idx.pbo.aw = ismember(T_this.aro_type, {'arousal', 'arousalemg'}) & ismember(T_this.ses, {'placebo'}) & T_this.is_awakening == 'true';
    pcfg.idx.pbo.cs = ismember(T_this.aro_type, {'arousal', 'arousalemg'}) & ismember(T_this.ses, {'placebo'}) & T_this.is_awakening == 'false';
    pcfg.idx.etc.aw = ismember(T_this.aro_type, {'arousal', 'arousalemg'}) & ismember(T_this.ses, {'etc120'}) & T_this.is_awakening == 'true';
    pcfg.idx.etc.cs = ismember(T_this.aro_type, {'arousal', 'arousalemg'}) & ismember(T_this.ses, {'etc120'}) & T_this.is_awakening == 'false';
    pcfg.idx.aro.aw = 695;
    pcfg.idx.aro.cs = 352;

    % Plot config
    pcfg.const = const;
    pcfg.xlim = [-50 20];
    pcfg.ylim_eeg =[-150 150];
    pcfg.ylim_sigma = [-3.25, 6.1];
    pcfg.ylim_pwr = [-0.15 0.67];
    pcfg.ylim_phase = [0 60];
    pcfg.margin = [0.06, 0.0, -0.035, -0.0];
    pcfg.resolution = 30;
    pcfg.nbins = 360/pcfg.resolution;
    pcfg.bins.edges = linspace(-pi, pi, pcfg.nbins+1);
    pcfg.bins.edges = pcfg.bins.edges - mean(diff(pcfg.bins.edges))./2;
    pcfg.bins.centers = pcfg.bins.edges(1:end-1) + diff(pcfg.bins.edges)/2;
    pcfg.bins.width = mean(diff(pcfg.bins.edges));
    pcfg.sigma_perm = Clust_this.y;

    close all
    % Circ stats
    Circ = run_circular_stats_tests(T_this, pcfg, roi, delay, Circ);

    % -------------------------------------------------------------------------
    % Create, assemble, and export Figure 2 for this region of interest
    plot_fig2_full(roi, T_this, Perms_this, Chans_this, chans, EEG, SIG, FSIG, pcfg);

end

disp('done')

%%
% Supplementary figure: phase-angle polar histograms at increasing delays
plot_supp_phase_coupling(T_this, pcfg);

%%
clear ans i pcfg

end
