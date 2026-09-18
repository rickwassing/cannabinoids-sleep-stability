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
% Init
H = struct();
cnt = 0;
FName = [];
XData = [];
YData = [];
% -------------------------------------------------------------------------
% Load and process files
for i = 1:length(SigmaFiles)
    % ---------------------------------------------------------------------
    % Load file and select parietal channels only
    SIGMA = LoadDataset(fullfile(SigmaFiles(i).folder, SigmaFiles(i).name), 'all');
    SIGMA = pop_select(SIGMA, 'channel', chans.idx);
    kv = filename2struct(SIGMA.setname);
    % ---------------------------------------------------------------------
    % Normalization factor
    normfact = median(SIGMA.data, 2, 'omitnan');
    % ---------------------------------------------------------------------
    % Extract bout indexes
    bouts = ([SIGMA.event(strcmpi({SIGMA.event.type}, 'boundary')).latency]);
    bouts = [0.5, bouts, SIGMA.pnts+0.5]; %#ok<AGROW>
    % For each bout
    for b = 1:length(bouts)-1
        % -----------------------------------------------------------------
        % Cut and append the data
        sig = pop_select(SIGMA, 'point', [bouts(b), bouts(b+1)]);
        sig.data = 10.*log10(sig.data./repmat(normfact, 1, sig.pnts));
        fsig = sig;
        fsig.data = detrend(fsig.data', 0, 'omitnan')'; % demean
        fsig = signalappend(fsig, 'zeros', 2, filtcfg);
        % -----------------------------------------------------------------
        % Prepend and append data to cover for filter edge artefact
        fsig = executeappending(fsig);
        % -----------------------------------------------------------------
        % Filter the data
        fsig = pop_firws(fsig, ...
            'fcutoff', filtcfg.cutoff, ...
            'ftype', 'bandpass', ...
            'wtype', filtcfg.wintype, ...
            'warg', 2, ...
            'forder', filtcfg.order, ...
            'plotfresp', false, ...
            'minphase', 0);
        % -----------------------------------------------------------------
        % Apply Hilbert
        fsig.cmx = hilbert(detrend(fsig.data', 0))';
        fsig.ang = angle(fsig.cmx);
        fsig.amp = abs(fsig.cmx);
        % -----------------------------------------------------------------
        % Remove appending
        fsig = executeappending(fsig, 'remove');
        % -----------------------------------------------------------------
        % Adjust time vector where zero indicates REM onset
        if b == 1
            remeplat = sig.event(find(strcmpi({sig.event.type}, 'remeps'), 1, 'first')).origlatency;
        else
            thislat = sig.event(find(strcmpi({sig.event.type}, 'remeps'), 1, 'first')).origlatency;
            remeplat = thislat - (remeplat+remepdur);
        end
        tzero = sig.event(find(strcmpi({sig.event.type}, 'remeps'), 1, 'first')).latency;
        remepdur = sig.event(find(strcmpi({sig.event.type}, 'remeps'), 1, 'first')).duration;
        sig.times = ((0:sig.pnts-1)-tzero)./sig.srate;
        % -----------------------------------------------------------------
        % Store
        cnt = cnt+1;
        H(cnt).sub = kv.sub;
        H(cnt).cond = kv.ses;
        H(cnt).rawsigma = sig.data';
        H(cnt).filtsigma = fsig.data';
        H(cnt).times = sig.times;
        H(cnt).cmx = fsig.cmx';
        H(cnt).ang = fsig.ang';
        H(cnt).amp = fsig.amp';
        H(cnt).remepdur = remepdur;
        H(cnt).remeplat = remeplat;
        % -----------------------------------------------------------------
        % Find peaks and troughs at the REM transition, and the ones before that
        H(cnt).event = prerempeakstroughs(H(cnt), evcfg);
    end
end
disp('done')

%% For each trace extract the peak amplitudes and duration

avsigma = [];
desmat = [];
for i = 1:length(H)

    tmp = table();
    tmp.sub = {H(i).sub};
    tmp.cond = {H(i).cond};
    tmp.s = {mean(H(i).rawsigma(1:3601, :), 2)};
    avsigma = [avsigma; tmp]; %#ok<AGROW>

    events = H(i).event;

    uuids = ascolumn(sort(arrayfun(@(x) getuuid(), 1:50, 'UniformOutput', false)));
    [eventGroups, eventKeys] = findgroups(events.id);

    tmp = table();
    tmp.latency = splitapply(@mean, events.latency+events.duration/2, eventGroups);
    tmp.duration = splitapply(@mean, events.duration/10, eventGroups);
    tmp.amplitude = splitapply(@mean, double(events.amplitude), eventGroups);
    tmp.stage = splitapply(@(x) x(1), events.stage, eventGroups);
    tmp.id = splitapply(@(x) x(1), events.id, eventGroups);
    tmp.sub = repmat({H(i).sub}, size(tmp, 1), 1);
    tmp.cond = repmat({H(i).cond}, size(tmp, 1), 1);
    tmp.remepdur = repmat(H(i).remepdur, size(tmp, 1), 1)./(60*10);
    tmp.remeplat = repmat(H(i).remeplat, size(tmp, 1), 1)./(60*10);
    tmp.id = uuids(findgroups(tmp.remepdur)); 

    desmat = [desmat; tmp]; %#ok<AGROW>
end

desmat.dur_z = zscore(desmat.duration);
desmat.amp_z = zscore(desmat.amplitude);
desmat.remdur_z = zscore(desmat.remepdur);
desmat.remlat_z = zscore(desmat.remeplat);

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

clear pcfg
pcfg.markersize = 5;

clear ids


close all;
Fig = figure('Color', 'w');
Fig.Units = 'centimeters';
Fig.Position = [1 12 18 8].*1;

clear Ax;
i = 0;

i = i+1;
Ax(i) = axes(...
    'NextPlot', 'add', ...
    'FontSize', 8, ...
    'Box', 'on', ...
    'LineWidth', 0.25, ...
    'Clipping', 'off', ...
    'Position', [0.05 0.77 0.4 0.18]);
plotHypnogram(Ax(i), SIG, 'LineWidth', 5, 'Hyp', HYP);

for b = 1:size(bouts, 1)
    YData = [-3.75; -3.75; 1.75; 1.75; -3.75];
    XData = [bouts(b, 1); bouts(b, 2); bouts(b, 2); bouts(b, 1); bouts(b, 1)];
    Vertices = [XData, YData];
    Faces = 1:4;
    patch(Ax(i), 'Faces', Faces, 'Vertices', Vertices, ...
        'EdgeColor', [0, 0, 0], ...
        'FaceColor', css_standard_colors('blue'), ...
        'FaceAlpha', 0.15, ...
        'EdgeAlpha', 0.30, ...
        'LineStyle', ':');
end

% Zoom lines
selbout = 2;
plot(Ax(i), [bouts(selbout, 1), 0], [-4.4 -6], ':', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5)
plot(Ax(i), [bouts(selbout, 2), Ax(i).XLim(2)], [-4.4 -6], ':', 'Color', [0.5, 0.5, 0.5], 'LineWidth', 0.5)

Ax(i).XTickLabel = {};
Ax(i).FontSize = 8;
Ax(i).XColor = 'w';
Ax(i).TickLength = [0, 0];

i = i + 1;
Ax(i) = plot_fig3_sigmatrace(Fig, H, pcfg);
Ax(i).Position = [0.05 0.55 0.4 0.18];

i = i + 1;
Ax(i) = plot_fig3_ampdur(desmat, 'duration');
Ax(i).Position = [0.55 0.6 0.13 0.29];

i = i + 1;
Ax(i) = plot_fig3_ampdur(desmat, 'amplitude');
Ax(i).Position = [0.82 0.6 0.13 0.29];

% Topoplot of selected channels
i = i+1;
% Create axes
Ax(i) = axes(Fig, 'NextPlot', 'add', 'Position', [0 0.175 0.1 0.175]);
topoplot(chans.idx, chans.locs, ...
    'style', 'blank', ...
    'electrodes', 'off', ...
    'emarker', {'.', 'k', 12, 1}, ...
    'emarkercolors', {[0, 0, 0]}, ...
    'hlinewidth', 1, ...
    'hcolor', [0.5, 0.5, 0.5], ...
    'colormap', Ax(2).UserData.CMap, ...
    'whitebk', 'on');

i = i + 1;
[Ax(i), leg] = plot_fig3_avsigmatrace(Fig, H, pcfg);
Ax(i).Position = [0.175 0.16 0.22 0.315];
leg.Position(1:2) = [0, Ax(i).Position(2)+Ax(i).Position(4)-leg.Position(4)];

i = i+1;
Ax(i) = plot_fig3_mdl(Fig, desmat, 'duration');
Ax(i).Position = [0.55 0.16 0.13 0.225];

i = i+1;
Ax(i) = plot_fig3_mdlhists(desmat, 'remeplat', 0:30:420, 'h');
Ax(i).Position = [sum(Ax(i-1).Position([1 3])), Ax(i-1).Position(2), 0.04, Ax(i-1).Position(4)];

i = i+1;
Ax(i) = plot_fig3_mdlhists(desmat, 'duration', 5:5:65, 'v');
Ax(i).Position = [Ax(i-2).Position(1), sum(Ax(i-2).Position([2 4])), Ax(i-2).Position(3), 0.09];

i = i+1;
Ax(i) = plot_fig3_mdl(Fig, desmat, 'amplitude');
Ax(i).Position = [0.82 0.16 0.13 0.225];

i = i+1;
Ax(i) = plot_fig3_mdlhists(desmat, 'remeplat', 0:30:420, 'h');
Ax(i).Position = [sum(Ax(i-1).Position([1 3])), Ax(i-1).Position(2), 0.04, Ax(i-1).Position(4)];

i = i+1;
Ax(i) = plot_fig3_mdlhists(desmat, 'amplitude', 0:0.2:3, 'v');
Ax(i).Position = [Ax(i-2).Position(1), sum(Ax(i-2).Position([2 4])), Ax(i-2).Position(3), 0.09];

plot_fig3_panellabels(Fig);

if strcmpi(chans.sel, 'fz')
    exportgraphics(Fig, './figures/figure3_fz.png', 'Resolution', 1200)
else
    exportgraphics(Fig, './figures/figure3_pz.png', 'Resolution', 1200)
end

