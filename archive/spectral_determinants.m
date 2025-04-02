% Load pre-arousal spectrograms and compare between continued sleep and
% state shifts
bidsdir = 'S:/Sleep/3. ACTIVE STUDIES/CUPID/Arousal paper backup_CANSLEEP';
COND = readtable('Scoring log and notes_CANSLEEP arousal.xlsx');
Taro = loadarousalcsvs(COND, 'selectcondition', false);

%%
min_interval = 15;
sleep_episodes = [-1, -2, 1, 2];
arousal_type = {'arousalemg'};

OUT = loadarousalepochs(Taro, min_interval, sleep_episodes, arousal_type, true);
disp('done loading')

%%

CMap = load('colormap_roma');
CMap = CMap.roma;
CLim = [-100, 150];
doLog = false;
doDetrend = false;

idx_ev = OUT.INDEP.episode == 2;
idx_ss = OUT.INDEP.stateshift == 1;

idx_time = OUT.spectimes >= -60 & OUT.spectimes <= 5;
idx_freq = OUT.specfreqs >= 0 & OUT.specfreqs <= 24;
%idx_freq = OUT.specfreqs >= 11 & OUT.specfreqs <= 16;

close all

Fig = figure('Position', [1000 761 300 180]);

XData = OUT.spectimes(idx_time);
YData = OUT.specfreqs(idx_freq);
CData = (squeeze(mean(mean(OUT.SPECDATA(:, idx_freq, idx_time, idx_ev & ~idx_ss), 4, 'omitnan'), 1, 'omitnan')));
if doLog
    if any(CData(:) < 0)
        CData = CData - min(CData(:));
    end
    CData = log10(CData);
end
if doDetrend
    mu = repmat(mean(CData(:, XData < -1), 2, 'omitnan'), 1, length(XData));
    CData =CData - mu;
end
Ax = plotarousalersp(Fig, XData, YData, CData, [XData(1), XData(end)], CMap, CLim);

exportgraphics(Fig, 'group-level/spectro_rem_fooof_long_contsleep.png', 'Resolution', 300)

Fig = figure('Position', [1300 761 300 180]);

XData = OUT.spectimes(idx_time);
YData = OUT.specfreqs(idx_freq);
CData = (squeeze(mean(mean(OUT.SPECDATA(:, idx_freq, idx_time, idx_ev & idx_ss), 4, 'omitnan'), 1, 'omitnan')));
if doLog
    if any(CData(:) < 0)
        CData = CData - min(CData(:));
    end
    CData = log10(CData);
end
if doDetrend
    mu = repmat(mean(CData(:, XData < -1), 2, 'omitnan'), 1, length(XData));
    CData =CData - mu;
end

Ax = plotarousalersp(Fig, XData, YData, CData, [XData(1), XData(end)], CMap, CLim);

exportgraphics(Fig, 'group-level/spectror_rem_fooof_long_stateshift.png', 'Resolution', 300)

%%

RES = struct();
RES.sws.beta = zeros(length(OUT.chanlocs), length(OUT.specfreqs));
RES.sws.tstat = zeros(length(OUT.chanlocs), length(OUT.specfreqs));
RES.sws.pval = zeros(length(OUT.chanlocs), length(OUT.specfreqs));
RES.asc.beta = zeros(length(OUT.chanlocs), length(OUT.specfreqs));
RES.asc.tstat = zeros(length(OUT.chanlocs), length(OUT.specfreqs));
RES.asc.pval = zeros(length(OUT.chanlocs), length(OUT.specfreqs));
RES.nrem.beta = zeros(length(OUT.chanlocs), length(OUT.specfreqs));
RES.nrem.tstat = zeros(length(OUT.chanlocs), length(OUT.specfreqs));
RES.nrem.pval = zeros(length(OUT.chanlocs), length(OUT.specfreqs));
RES.rem.beta = zeros(length(OUT.chanlocs), length(OUT.specfreqs));
RES.rem.tstat = zeros(length(OUT.chanlocs), length(OUT.specfreqs));
RES.rem.pval = zeros(length(OUT.chanlocs), length(OUT.specfreqs));

t = now;

idx_time = OUT.spectimes >= -10 & OUT.spectimes <= -1;
for episode = [2 3]
    idx_ev = OUT.INDEP.episode == episode;
    switch episode
        case -2
            fld = 'sws';
        case 1
            fld = 'asc';
        case 2
            fld = 'rem';
        case 3
            fld = 'nrem';
            idx_ev = OUT.INDEP.episode == -2 |  OUT.INDEP.episode == 1;
    end
    for idx_chan = 1:length(OUT.chanlocs)
        for idx_freq = 1:length(OUT.specfreqs)-1 % FOOOF does not have first freq
            
            T = table();
            T.Y = double(ascolumn(squeeze(mean(OUT.SPECDATA(idx_chan, idx_freq, idx_time, idx_ev), 3, 'omitnan'))));
            T.X = OUT.INDEP.stateshift(idx_ev);
            T.sub = OUT.INDEP.subject(idx_ev);

            mdl = fitglme(T, 'Y ~ 1 + X + (1|sub)', 'DummyVarCoding', 'effects');
           
            RES.(fld).beta(idx_chan, idx_freq) = mdl.Coefficients.Estimate(2);
            RES.(fld).tstat(idx_chan, idx_freq) = mdl.Coefficients.tStat(2);
            RES.(fld).pval(idx_chan, idx_freq) = mdl.Coefficients.pValue(2);
    
            t = remainingTime(t, 2*length(OUT.chanlocs)*(length(OUT.specfreqs)-1));
        end
    end
end

%%

close all

idx_time = OUT.spectimes >= -15 & OUT.spectimes <= -1;
idx_freq = OUT.specfreqs > 0 & OUT.specfreqs <= 24;

XData = OUT.specfreqs(idx_freq);

doLog = true;

YLim = [1.75, 3.2];

for episode = 2%[-2, 1, 2]

    idx_ev = OUT.INDEP.episode == episode;
    idx_ss = OUT.INDEP.stateshift == 1;

    switch episode
        case -2
            fld = 'sws';
            highlight = [];
        case 1
            fld = 'asc';
            highlight = [...
                13, 15.5; ...
                ];
        case 2
            fld = 'rem';
            highlight = [];
        case 3
            fld = 'nrem';
            idx_ev = OUT.INDEP.episode < 2;
    end

    CData_A = mean(mean(OUT.SPECDATA(:, :, idx_time, idx_ev & ~idx_ss), 4, 'omitnan'), 3, 'omitnan');
    CData_B = mean(mean(OUT.SPECDATA(:, :, idx_time, idx_ev & idx_ss), 4, 'omitnan'), 3, 'omitnan');

    if doLog
        CData_A = log10(CData_A')';
        CData_B = log10(CData_B')';
    end
    Fig = plot2dspectrum(XData, CData_A, CData_B, RES.(fld).pval(:, idx_freq), CMap, YLim, idx_ev, idx_ss, highlight);
    
    exportgraphics(Fig, ['group-level/avspectro_', fld, '.png'], 'Resolution', 600)
end

%% Phase angle of arousals

idx_ev = ...
    (Taro.sleepepisode == -2 | Taro.sleepepisode == 1) & ...
    Taro.interval >= 15 & ...
    strcmpi(Taro.newtype, 'arousalemg');

phi = Taro.ism_phs(idx_ev);
idx_ss = Taro.stateshift(idx_ev) == 1;
idx_ss(isnan(phi)) = [];
phi(isnan(phi)) = [];

pval = circ_wwtest(phi, idx_ss);

close all

Fig = figure('Position', [896 683 566 223]);

Ax = axes(...
    'Box', 'on', ...
    'NextPlot', 'add', ...
    'Color', [233, 236, 239]./255, ...
    'XColor', [0.5, 0.5, 0.5], ...
    'YColor', [0.5, 0.5, 0.5], ...
    'XGrid', 'on', ...
    'XTick', pi/2:pi/2:2*pi+pi/2, ...
    'XTickLabel', {'90° (falling)', '180° (trough)', '270° (rising)', '0° (peak)'}, ...
    'YTick', [-.9, .9], ...
    'YTickLabel', {'low sigma power', 'high sigma power'}, ...
    'YGrid', 'off', ...
    'YMinorGrid', 'off', ...
    'GridAlpha', 0.25, ...
    'LineWidth', 0.5, ...
    'TickLength', [0, 0], ...
    'FontSize', 10, ...
    'Layer', 'top', ...
    'Units', 'normalized', ...
    'Toolbar', [], ...
    'Interactions', []);

[mu, ul, ll] = circ_mean(phi(~idx_ss));
XData = phi(~idx_ss);
XData = [XData; XData + 2.*pi];
YData = cos(XData) + 0.1;

clear h

scatter(XData, YData, 'o', ...
    'MarkerEdgeColor', CMap(1, :), ...
    'MarkerFaceColor', CMap(1, :), ...
    'MarkerEdgeAlpha', 0.05, ...
    'MarkerFaceAlpha', 0.05);

plot((ll:pi/50:ul)+2*pi, cos(ll:pi/50:ul)+0.1, '-w', 'LineWidth', 3)

h(1) = plot(mu+2*pi, cos(mu) + 0.1, 'o', ...
    'MarkerFaceColor', CMap(1, :), ...
    'MarkerEdgeColor', 'w', ...
    'MarkerSize', 6);

[mu, ul, ll] = circ_mean(phi(idx_ss));
XData = phi(idx_ss);
XData = [XData; XData + 2.*pi];
YData = cos(XData) - 0.1;

scatter(XData, YData, 'o', ...
    'MarkerEdgeColor', CMap(end, :), ...
    'MarkerFaceColor', CMap(end, :), ...
    'MarkerEdgeAlpha', 0.05, ...
    'MarkerFaceAlpha', 0.05);


plot((ll:pi/50:ul)+2*pi, cos(ll:pi/50:ul)-0.1, '-w', 'LineWidth', 3)

h(2) = plot(mu+2*pi, cos(mu) - 0.1, 'o', ...
    'MarkerFaceColor',  CMap(end, :), ...
    'MarkerEdgeColor', 'w', ...
    'MarkerSize', 6);

Ax.XLim = [pi/2, 2.*pi+pi/2];
Ax.YLim = [-1.2 1.2];

leg = legend(h, {sprintf('Continued sleep (N = %i)', sum(~idx_ss)), sprintf('State-shift (N = %i)', sum(idx_ss))}, 'Box', 'off', 'FontSize', 10);
leg.Location = 'southeast';

%% PLOT HR versus SIGMA POW

sigma_filepath = 'derivatives/EEG-preproc/sub-r0091/ses-1/sub-r0091_ses-1_task-psg_run-1_desc-sigma_eeg.set';
fit_filepath = 'derivatives/EEG-preproc/sub-r0091/ses-1/sub-r0091_ses-1_task-psg_run-1_desc-ismfit_struct.mat';
SIGMA = LoadDataset(sigma_filepath, 'all'); 
hr = SIGMA.data(end, :);
filtorder = pop_firwsord('hamming', SIGMA.srate, 0.01);
INFRA = pop_firws(SIGMA, ...
    'fcutoff', [0.005, 0.055], ...
    'ftype', 'bandpass', ...
    'wtype', 'hamming', ...
    'forder', filtorder, ...
    'minphase', 0);
idx_nan = sum(isnan(INFRA.data), 1) > 0;
INFRA.data(:, idx_nan) = 0;
FIT = load(fit_filepath);
FIT = FIT.FIT;
idx_chan = [FIT.peak] > 1;
idx_chan(end) = []; % remove the last element (the average)
infra_mu = mean(INFRA.data(idx_chan, :), 1);
tmp = hilbert(infra_mu);
amp = abs(tmp)';
phs = angle(tmp)';
infra_mu(idx_nan) = NaN;
amp(idx_nan) = NaN;
phs(idx_nan) = NaN;

%%

Fig = figure('Color', 'w')';
hold on
idx_chan = 21;
plot(SIGMA.times, nanzscore(mean(SIGMA.data(idx_chan, :))))
plot(SIGMA.times, nanzscore(infra_mu)+2, '-', 'LineWidth', 1.5)

for i = 1:length(SIGMA.event)
    if ~strcmpi(SIGMA.event(i).type, 'arousalemg')
        continue
    end
    XData = [SIGMA.event(i).latency, SIGMA.event(i).latency+SIGMA.event(i).duration] ./ SIGMA.srate;
    YData = [7, 7];
    plot(XData, YData, '-', 'Color', CMap(end, :), 'LineWidth', 4)
end
%plot(SIGMA.times, zscore(INFRA.data(end, :)))


%% PLOT SIGMA RESPONSE

idx_ev = OUT.INDEP.episode < 2 & OUT.INDEP.interval >= 1;
idx_ss = OUT.INDEP.stateshift == 1;
idx_time = OUT.spectimes >= -60 & OUT.spectimes <= 15;
idx_freq = OUT.specfreqs >= 11 & OUT.specfreqs <= 16;

figure

clear Ax h

Ax(1) = axes(...
    'Box', 'on', ...
    'NextPlot', 'add', ...
    'Color', [233, 236, 239]./255, ...
    'XColor', [0.5, 0.5, 0.5], ...
    'YColor', [0.5, 0.5, 0.5], ...
    'XGrid', 'on', ...
    'YGrid', 'off', ...
    'YMinorGrid', 'off', ...
    'GridAlpha', 0.25, ...
    'LineWidth', 0.5, ...
    'TickLength', [0, 0], ...
    'FontSize', 10, ...
    'Layer', 'top', ...
    'Units', 'normalized', ...
    'Toolbar', [], ...
    'Interactions', []);

XData = OUT.spectimes(idx_time);

YData = squeeze(mean(mean(mean(OUT.SPECDATA(:, idx_freq, idx_time, idx_ev & ~idx_ss), 4, 'omitnan'), 2, 'omitnan'), 1, 'omitnan'));
tstat = tinv(0.975, sum(idx_ev & ~idx_ss));
SData = std(squeeze(mean(mean(OUT.SPECDATA(:, idx_freq, idx_time, idx_ev & ~idx_ss), 1, 'omitnan'), 2, 'omitnan')), [], 2, 'omitnan');
SData = tstat .* SData ./ sum(idx_ev & ~idx_ss);

patch('XData', [asrow(XData), fliplr(asrow(XData))], 'YData', [asrow(YData), fliplr(asrow(YData))] + [asrow(SData), -1.*fliplr(asrow(SData))], ...
    'LineStyle', 'none', ...
    'FaceColor', CMap(1, :), ...
    'FaceAlpha', 0.5)
h(1) = plot(XData, YData, '-', 'Color', CMap(1, :));

YData = squeeze(mean(mean(mean(OUT.SPECDATA(:, idx_freq, idx_time, idx_ev & idx_ss), 4, 'omitnan'), 2, 'omitnan'), 1, 'omitnan'));
tstat = tinv(0.975, sum(idx_ev & idx_ss));
SData = std(squeeze(mean(mean(OUT.SPECDATA(:, idx_freq, idx_time, idx_ev & idx_ss), 1, 'omitnan'), 2, 'omitnan')), [], 2, 'omitnan');
SData = tstat .* SData ./ sum(idx_ev & idx_ss);

patch('XData', [asrow(XData), fliplr(asrow(XData))], 'YData', [asrow(YData), fliplr(asrow(YData))] + [asrow(SData), -1.*fliplr(asrow(SData))], ...
    'LineStyle', 'none', ...
    'FaceColor', CMap(end, :), ...
    'FaceAlpha', 0.5)
h(2) = plot(XData, YData, '-', 'Color', CMap(end, :));

leg = legend(h, {sprintf('Continued sleep (N = %i)', sum(idx_ev & ~idx_ss)), sprintf('State-shift (N = %i)', sum(idx_ev & idx_ss))}, 'Box', 'off', 'FontSize', 10);
leg.Location = 'northwest';

Ax.XLim = [min(OUT.spectimes(idx_time)), max(OUT.spectimes(idx_time))];
Ax.XLabel.String = 'time to arousal (s)';
Ax.YLabel.String = 'Sigma power (a.u.)';


%% PLOT HR RESPONSE

idx_ev = OUT.INDEP.episode == 2 & OUT.INDEP.interval >= 1;
idx_ss = OUT.INDEP.stateshift == 1;

Fig = figure();

clear Ax

Ax(1) = axes(...
    'Box', 'on', ...
    'NextPlot', 'add', ...
    'Color', [233, 236, 239]./255, ...
    'XColor', [0.5, 0.5, 0.5], ...
    'YColor', [0.5, 0.5, 0.5], ...
    'XGrid', 'on', ...
    'YGrid', 'off', ...
    'YMinorGrid', 'off', ...
    'GridAlpha', 0.25, ...
    'LineWidth', 0.5, ...
    'TickLength', [0, 0], ...
    'FontSize', 10, ...
    'Layer', 'top', ...
    'Units', 'normalized', ...
    'Toolbar', [], ...
    'Interactions', []);

XData = -60:0.2:14.8;

YData = mean(OUT.HR(:, idx_ev & ~idx_ss), 2, 'omitnan');
tstat = tinv(0.975, sum(idx_ev & ~idx_ss));
SData = tstat .* (std(OUT.HR(:, idx_ev & ~idx_ss), [], 2, 'omitnan') ./ sum(idx_ev & ~idx_ss));

patch('XData', [XData, fliplr(XData)], 'YData', [asrow(YData), fliplr(asrow(YData))] + [asrow(SData), -1.*fliplr(asrow(SData))], ...
    'LineStyle', 'none', ...
    'FaceColor', CMap(1, :), ...
    'FaceAlpha', 0.5)
h(1) = plot(XData, YData, '-', 'Color', CMap(1, :));

YData = mean(OUT.HR(:, idx_ev & idx_ss), 2, 'omitnan');
tstat = tinv(0.975, sum(idx_ev & idx_ss));
SData = tstat .* (std(OUT.HR(:, idx_ev & idx_ss), [], 2, 'omitnan') ./ sum(idx_ev & idx_ss));

patch('XData', [XData, fliplr(XData)], 'YData', [asrow(YData), fliplr(asrow(YData))] + [asrow(SData), -1.*fliplr(asrow(SData))], ...
    'LineStyle', 'none', ...
    'FaceColor', CMap(end, :), ...
    'FaceAlpha', 0.5)
h(2) = plot(XData, YData, '-', 'Color', CMap(end, :));

leg = legend(h, {sprintf('Continued sleep (N = %i)', sum(idx_ev & ~idx_ss)), sprintf('State-shift (N = %i)', sum(idx_ev & idx_ss))}, 'Box', 'off', 'FontSize', 10);
leg.Location = 'northwest';

Ax.XLim = [-60, 5];
Ax.XLabel.String = 'time to arousal (s)';
Ax.YLabel.String = 'Heart rate (BPM)';

exportgraphics(Fig, 'group-level/heartrate-response_rem.png', 'Resolution', 600)

%%

labels = {'sigma'};
bands = [...
    13, 15.5; ...
    ];

T = now();

for i = 1:size(Taro, 1)

    T = remainingTime(T, size(Taro, 1), 'simple', true);
    
    kv = filename2struct(Taro.filename{i});
    kv = rmfield(kv, 'run');
    kv = rmfield(kv, 'ses');
    kv.desc = sprintf('e%i', Taro.id(i));
    kv.filetype = 'powerspect.mat';

    psd = struct();
    psd.filename = struct2filename(kv);
    psd.filepath = [bidsdir, '/derivatives/EEG-output-fstlvl/sub-', kv.sub];
    if exist(fullfile(psd.filepath, psd.filename), 'file') == 2
        continue
    end
    psd.subject = kv.sub;
    psd.session = '';
    psd.task = kv.task;
    psd.group = '';
    psd.condition = ifelse(Taro.stateshift(i) == 1, 'stateshift', 'cont-sleep');
    psd.nbchan = 178;
    psd.trials = 1;
    psd.data = squeeze(TF.psddata(:, :, i));
    psd.freqs = TF.psdfreqs;
    psd.freqstep = mean(diff(TF.psdfreqs));
    psd.bands = struct();
    for j = 1:size(bands, 1)
        idx_f = psd.freqs >= bands(j, 1) & psd.freqs <= bands(j, 2);
        psd.bands(j).label = labels{j};
        psd.bands(j).type = 'absolute';
        psd.bands(j).freqrange = bands(j, :);
        psd.bands(j).data = squeeze(mean(TF.psddata(:, idx_f, i), 2, 'omitnan'));
    end
    psd.chanlocs = readlocs('GSN-HydroCel-257.sfp');
    psd.chanlocs(~ismember({psd.chanlocs.labels}, TF.chans)) = [];
    psd.etc = struct();
    SaveDataset(psd, 'matrix');
end
app = struct();
app.State.Protocol.Path = pwd;
app.State.Subjects = GetSubjects(app);
app.State.Files = GetFiles(app);
