clear
clc

Files = dir('derivatives/aroepochs/sub-*/ses-*/sub*spindlesaroboutnrem2*.mat');

i=1;
EEG = load(fullfile(Files(i).folder, Files(i).name));
ARO = EEG.ARO;

AData = struct(); % Structure to save indices of arousals occuring pre-arousal
AData.cs = nan(0, length(ARO.times));
AData.aw = nan(0, length(ARO.times));

CData = struct(); % Structure to save the prearousal timeseries to
CData.cs = nan(size(ARO.data, 1), length(ARO.times(1):0.5:ARO.times(end)), 0, 'single');
CData.aw = nan(size(ARO.data, 1), length(ARO.times(1):0.5:ARO.times(end)), 0, 'single');

Subj = struct(); % To save the subject IDs
Subj.cs = [];
Subj.aw = [];

Cond = struct(); % To save the conditions (PLB vs ETC)
Cond.cs = [];
Cond.aw = [];

rt = now();
for i = 1:length(Files)
    kv = filename2struct(Files(i).name);
    % Load
    EEG = load(fullfile(Files(i).folder, Files(i).name));
    ARO = EEG.ARO;
    clear EEG;
    % -------------------------------------------------------------------------
    % Extract awakening and continued sleep trials
    idx_aw = ARO.event.is_awakening;
    AData.cs = cat(1, AData.cs, ARO.aro(~idx_aw, :));
    AData.aw = cat(1, AData.aw, ARO.aro(idx_aw, :));
    % -------------------------------------------------------------------------
    % Average across time (half seconds)
    d = nan(size(ARO.data, 1), length(ARO.times(1):0.5:ARO.times(end)), size(ARO.data, 3), 'single');
    cnt = 0;
    for t = ARO.times(1):0.5:ARO.times(end)
        idx_t = find(ARO.times >= (t-0.25) & ARO.times < t+0.25);
        cnt = cnt+1;
        d(:, cnt, :) = round(mean(ARO.data(:, idx_t, :), 2, 'omitnan'));
    end
    % -------------------------------------------------------------------------
    CData.cs = cat(3, CData.cs, d(:, :, ~idx_aw));
    CData.aw = cat(3, CData.aw, d(:, :, idx_aw));
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
foutname = sprintf('group-level/%s_prearo-spindlesnrem2.mat', datestr(now, 'yyyymmddTHHMM'));
save(foutname, 'AData', '-mat', '-v7.3');
save(foutname, 'CData', '-append');
save(foutname, 'Cond', '-append');
save(foutname, 'Subj', '-append');

%% APPLY GLM
% -------------------------------------------------------------------------
% Get channel location clusters
ARO.chanlocs = channel_clusters(ARO.chanlocs, 'MFF');
% -------------------------------------------------------------------------
% Init table for modelling
tbl = table();
tbl.is_awake = [true(size(CData.aw, 3), 1); false(size(CData.cs, 3), 1)];
tbl.cond = [Cond.aw; Cond.cs];
tbl.sub = [Subj.aw; Subj.cs];
% -------------------------------------------------------------------------
% Output structure
GLM = struct();
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Linear mixed model, IS_SPINDLE ~ 1 + IS_AWAKE*COND
GLM.spd.awake.tstat = nan(178, 180);
GLM.spd.awake.pval = nan(178, 180);
GLM.spd.cond.tstat = nan(178, 180);
GLM.spd.cond.pval = nan(178, 180);
GLM.spd.intx.tstat = nan(178, 180);
GLM.spd.intx.pval = nan(178, 180);
% -------------------------------------------------------------------------
% Output structure
PostHoc = struct();
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Linear mixed model, IS_SPINDLE ~ 1 + IS_AWAKE (separate for placebo and ETC120)
PostHoc.plc.awake.tstat = nan(178, 180);
PostHoc.plc.awake.pval = nan(178, 180);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
PostHoc.etc.awake.tstat = nan(178, 180);
PostHoc.etc.awake.pval = nan(178, 180);

% -------------------------------------------------------------------------
rt = now(); %#ok<TNOW1>
for c = 1:178
    % ---------------------------------------------------------------------
    % Get channel indices
    idx_c = c; % [ARO.chanlocs.cluster] == c;
    for t = 1:180
        % -------------------------------------------------------------
        % Replace power values in the table
        tbl.depvar = double([...
            squeeze(CData.aw(idx_c, t, :)); ...
            squeeze(CData.cs(idx_c, t, :))]);
        % -------------------------------------------------------------
        % Fit full linear mixed model
        try
            glme = fitglme(tbl, 'depvar ~ 1 + is_awake*cond + (1|sub)', 'DummyVarCoding', 'effects', 'Distribution', 'Binomial');
            % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            GLM.spd.awake.tstat(c, t) = glme.Coefficients.tStat(2);
            GLM.spd.awake.pval(c, t) = glme.Coefficients.pValue(2);
            GLM.spd.cond.tstat(c, t) = glme.Coefficients.tStat(3);
            GLM.spd.cond.pval(c, t) = glme.Coefficients.pValue(3);
            GLM.spd.intx.tstat(c, t) = glme.Coefficients.tStat(4);
            GLM.spd.intx.pval(c, t) = glme.Coefficients.pValue(4);
        catch
            GLM.spd.awake.tstat(c, t) = nan;
            GLM.spd.awake.pval(c, t) = nan;
            GLM.spd.cond.tstat(c, t) = nan;
            GLM.spd.cond.pval(c, t) = nan;
            GLM.spd.intx.tstat(c, t) = nan;
            GLM.spd.intx.pval(c, t) = nan;
        end
        % -------------------------------------------------------------
        % Fit linear mixed model for placebo
        try
            idx_rows = strcmpi(tbl.cond, 'placebo');
            glme_plc = fitglme(tbl(idx_rows, :), 'depvar ~ 1 + is_awake + (1|sub)', 'DummyVarCoding', 'effects', 'Distribution', 'Binomial');
            % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            PostHoc.plc.awake.tstat(c, t) = glme_plc.Coefficients.tStat(2);
            PostHoc.plc.awake.pval(c, t) = glme_plc.Coefficients.pValue(2);
        catch
            PostHoc.plc.awake.tstat(c, t) = nan;
            PostHoc.plc.awake.pval(c, t) = nan;
        end
        % -------------------------------------------------------------
        % Fit linear mixed model for etc120
        try
            idx_rows = strcmpi(tbl.cond, 'etc120');
            glme_etc = fitglme(tbl(idx_rows, :), 'depvar ~ 1 + is_awake + (1|sub)', 'DummyVarCoding', 'effects', 'Distribution', 'Binomial');
            % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
            PostHoc.etc.awake.tstat(c, t) = glme_etc.Coefficients.tStat(2);
            PostHoc.etc.awake.pval(c, t) = glme_etc.Coefficients.pValue(2);
        catch
            PostHoc.etc.awake.tstat(c, t) = nan;
            PostHoc.etc.awake.pval(c, t) = nan;
        end
        % -------------------------------------------------------------
        % Display remaining time
        rt = remainingTime(rt, 178*180, 'simple');
    end
end

% RUN THIS AGAIN TO SAVE 'POSTHOC'
save(sprintf('group-level/%s_modeldata-spindlesnrem2.mat', datestr(now, 'yyyymmddTHHMM')), 'GLM', 'PostHoc', 'glme', 'glme_plc', 'glme_etc', '-mat', '-v7.3') %#ok<DATST,TNOW1>

%% PRE-AROUSAL SPECTRUM
% -------------------------------------------------------------------------
% Parameters
times = -90:0.5:30;
times_orig = -90:1/256:30;
XLim = [-30, 10];
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Change here for the different conditions
fld_cond = 'etc';
switch fld_cond
    case 'plc'
        idx_cond_cs = strcmpi(Cond.cs, 'placebo');
        idx_cond_aw = strcmpi(Cond.aw, 'placebo');
    case 'etc'
        idx_cond_cs = ~strcmpi(Cond.cs, 'placebo');
        idx_cond_aw = ~strcmpi(Cond.aw, 'placebo');
end
% -------------------------------------------------------------------------
% Create figure
Fig = figure();
Fig.Units = 'centimeters';
Fig.Position = [5, 5, 9, 10];
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
    'XData', [times_orig(1), times_orig, times_orig(end)], ...
    'YData', [0, mean(AData.cs(idx_cond_cs, :), 1, 'omitnan'), 0], ...
    'LineStyle', 'none', ...
    'FaceColor', standard_colors('blue'), ...
    'FaceAlpha', 0.25)
patch(...
    'XData', [times_orig(1), times_orig, times_orig(end)], ...
    'YData', [0, mean(AData.aw(idx_cond_aw, :), 1, 'omitnan'), 0], ...
    'LineStyle', 'none', ...
    'FaceColor', standard_colors('red'), ...
    'FaceAlpha', 0.25)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
YTick = max([mean(AData.cs(idx_cond_cs, times_orig < 0), 1, 'omitnan'), mean(AData.aw(idx_cond_aw, times_orig < 0), 1, 'omitnan')]);
Ax.YTick = [ceil(YTick*100)/100, 1];

% -------------------------------------------------------------------------
% Pre-arousal spindles
Ax = axes(Fig);
Ax.NextPlot = 'add';
Ax.TickDir = 'out';
Ax.Layer = 'top';
Ax.Box = 'on';
Ax.XLim = XLim;
Ax.XTick = sort([XLim, 0]);
Ax.YLim = [0, 0.5];
Ax.XTick = sort([XLim, 0]);
Ax.YTick = [0, 1];
Ax.XTickLabel = {''};
Ax.YGrid = 'on';
Ax.Position = [0.3 0.3 0.55 0.4];
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Extract the pre-arousal spectra of the two state-shift conditions
AvCDataCS = CData.cs(:, :, idx_cond_cs);
AvCDataAW = CData.aw(:, :, idx_cond_aw);
AvCDataCS = mean(AvCDataCS, 1, 'omitnan'); % average across channels
AvCDataAW = mean(AvCDataAW, 1, 'omitnan'); % average across channels
AvCDataCS = mean(AvCDataCS, 3, 'omitnan'); % average across trials
AvCDataAW = mean(AvCDataAW, 3, 'omitnan'); % average across trials
plot(times, squeeze(AvCDataCS), '-', 'Color', standard_colors('blue'))
plot(times, squeeze(AvCDataAW), '-', 'Color', standard_colors('red'))
plot([0, 0], [0, 31], '-w')
% % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% XData = -90:0.5:-0.5;
% YData = 1:30;
% n_tests = 60*30; % number of correlated tests (60 timepoints and 30 freqs)
% alpha_sidak = (1-normcdf(3.1))*2; %1-(1-0.05).^(1/n_tests);
% alpha_posthoc = (1-normcdf(1.96)); %1-(1-0.05).^(1/n_tests);
% switch sleepstage
%     case 'nrem2'
%         MainSData = double(squeeze(sum(GLM.logpow.awake.pval < alpha_sidak, 1)) > 0)';
%         contour(XData, YData, MainSData, 1, 'k');
%         IntxSData = double(squeeze(sum(GLM.logpow.intx.pval < alpha_sidak, 1)) > 0)';
%         contour(XData, YData, IntxSData, 1, 'w', 'LineWidth', 1.5);
%         PostHocSData = double(squeeze(sum(PostHoc.(fld_cond).awake.pval < alpha_posthoc, 1)) > 0)';
%         PostHocSData = PostHocSData.*IntxSData;
%         contour(XData, YData, PostHocSData, 1, 'k');
%     case 'rem'
%         IntxSData = double(squeeze(sum(GLM.logpow.intx.pval < alpha_sidak, 1)) > 0)';
%         contour(XData, YData, IntxSData, 1, 'w', 'LineWidth', 1.5);
%         MainSData = double(squeeze(sum(GLM.rem.awake.pval < alpha_sidak, 1)) > 0)';
%         contour(XData, YData, MainSData, 1, 'k');
% end
% % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% switch sleepstage
%     case 'nrem2'
%         plot(Ax, [-12.75, -4.25, -4.25, -12.75, -12.75], [14.5, 14.5, 16.5, 16.5, 14.5], '-', 'Color', standard_colors('green'))
%         plot(Ax, [-8.75, -4.25, -4.25, -8.75, -8.75], [10.5, 10.5, 14.5, 14.5, 10.5], '-', 'Color', standard_colors('green'))
%     case 'rem'
%         plot(Ax, [-18.5, -11.25, -11.25, -18.5, -18.5], [15.5, 15.5, 24.5, 24.5, 15.5], '-', 'Color', standard_colors('green'))
%         plot(Ax, [-29, -9.75, -9.75, -29, -29], [11.5, 11.5, 15.5, 15.5, 11.5], '-', 'Color', standard_colors('green'))
%         plot(Ax, [-15.25, -8.5, -8.5, -15.25, -15.25], [6.5, 6.5, 11.5, 11.5, 6.5], '-', 'Color', standard_colors('green'))
% end
% % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Ax.Colormap = CMapRoma;

% -------------------------------------------------------------------------
% Number of significant channels
Ax = axes(Fig);
Ax.NextPlot = 'add';
Ax.TickDir = 'out';
Ax.Layer = 'top';
Ax.Box = 'on';
Ax.XLim = XLim;
Ax.XTick = sort([XLim, 0]);
Ax.YLim = [0, 1];
Ax.YTick = [0, 1];
Ax.XGrid = 'on';
Ax.Position = [0.3 0.15 0.55 0.15];
Ax.YAxisLocation = 'right';
Ax.FontSize = 8;
Ax.XLabel.String = 'Time (s)';
Ax.XLabel.FontSize = 8;
Ax.XLabel.Position(2) = 0.000045;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
XData = times(1:180);
YData = sum(GLM.spd.awake.pval <= 0.05) ./ size(GLM.spd.awake.pval, 1);
plot(Ax, XData, YData, ...
    'Color', 'k', ...
    'LineWidth', 2);
YData = sum(GLM.spd.cond.pval <= 0.05) ./ size(GLM.spd.cond.pval, 1);
plot(Ax, XData, YData, ...
    'Color', standard_colors('red'), ...
    'LineWidth', 2);
YData = sum(GLM.spd.intx.pval <= 0.05) ./ size(GLM.spd.intx.pval, 1);
plot(Ax, XData, YData, ...
    'Color', standard_colors('green'), ...
    'LineWidth', 2);

% 
% switch sleepstage
%     case 'nrem2'
%         figoutname = sprintf('figures/fig_2_%s_%s.png', sleepstage, fld_cond);
%     case 'rem'
%         figoutname = sprintf('figures/fig_2_%s.png', sleepstage);
% end
% exportgraphics(Fig, figoutname, 'Resolution', 1200)

