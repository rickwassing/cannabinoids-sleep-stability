function plot_fig2_full(roi, T_this, Perms_this, Chans_this, chans, EEG, SIG, FSIG, pcfg)
% -------------------------------------------------------------------------
% Assemble and export Figure 2 for a single region of interest ('fz' or
% 'pz'): EEG/sigma traces (fz only), channel topoplot, averaged prearousal
% sigma traces with Pr(arousal), phase-angle polar histograms, and
% instantaneous ISF amplitude.
% -------------------------------------------------------------------------
Fig = figure('Color', 'w');
Fig.Units = 'centimeters';
Fig.Position = [1 12 18 7].*1;
% -------------------------------------------------------------------------
% Init axes array
clear Ax
i = 0;
% -------------------------------------------------------------------------
% Plot EEG trace for panel A
if strcmpi(roi, 'fz')
    i = i+1;
    Ax(i) = plot_fig2_eegtrace(Fig, EEG, 'aw', pcfg); %#ok<*SAGROW>
    Ax(i).Position = [-0.015 1-0.2 0.47 0.2]+pcfg.margin;
    % Plot sigma power
    i = i+1;
    Ax(i) = plot_fig2_sigmatrace(Fig, SIG, FSIG, 'aw', pcfg);
    Ax(i).Position = [-0.015 1-0.37 0.47 0.175]+pcfg.margin;
end
% -------------------------------------------------------------------------
% Plot EEG trace for panel B
if strcmpi(roi, 'fz')
    i = i+1;
    Ax(i) = plot_fig2_eegtrace(Fig, EEG, 'cs', pcfg);
    Ax(i).Position = [0.48 1-0.2 0.47 0.2]+pcfg.margin;
    % Plot sigma power
    i = i+1;
    Ax(i) = plot_fig2_sigmatrace(Fig, SIG, FSIG, 'cs', pcfg);
    Ax(i).Position = [0.48 1-0.37 0.47 0.175]+pcfg.margin;
end
% -------------------------------------------------------------------------
% Topoplot of selected channels
i = i+1;
% Create axes
Ax(i) = axes(Fig, 'NextPlot', 'add', 'Position', [0 0.175 0.1 0.225]);
topoplot(Chans_this, chans.locs, ...
    'style', 'blank', ...
    'electrodes', 'off', ...
    'emarker', {'.', 'k', 12, 1}, ...
    'emarkercolors', {[0, 0, 0]}, ...
    'hlinewidth', 1, ...
    'hcolor', [0.5, 0.5, 0.5], ...
    'whitebk', 'on');

% -------------------------------------------------------------------------
% Plot averaged prearousal sigma timeseries between CS and AW arousals for
% AWAKENINGS AROUSALS
i = i+1;
Ax(i) = plot_fig2_avsigmatrace(Fig, T_this, Perms_this, 'aw', pcfg);
Ax(i).Position = [0.11 0.22 0.18 0.2]+pcfg.margin;
% -------------------------------------------------------------------------
% Plot Pr Arousals
i = i+1;
Ax(i) = plot_fig2_probarousal(Fig, T_this, Perms_this, 'aw', pcfg);
Ax(i).Position = [Ax(i-1).Position(1), Ax(i-1).Position(2)+Ax(i-1).Position(4), Ax(i-1).Position(3) 0.1];
Ax(i).XLim = Ax(i-1).XLim;
Ax(i).XTick = Ax(i-1).XTick;

% -------------------------------------------------------------------------
% Plot averaged prearousal sigma timeseries between CS and AW arousals for
% CONTINUED SLEEP AROUSALS
i = i+1;
Ax(i) = plot_fig2_avsigmatrace(Fig, T_this, Perms_this, 'cs', pcfg);
Ax(i).Position = [0.28 0.22 0.18 0.2]+pcfg.margin;
% -------------------------------------------------------------------------
% Plot Pr Arousals
i = i+1;
Ax(i) = plot_fig2_probarousal(Fig, T_this, Perms_this, 'cs', pcfg);
Ax(i).Position = [Ax(i-1).Position(1), Ax(i-1).Position(2)+Ax(i-1).Position(4), Ax(i-1).Position(3) 0.1];
Ax(i).XLim = Ax(i-1).XLim;
Ax(i).XTick = Ax(i-1).XTick;

% -------------------------------------------------------------------------
% Phase angle of AWAKENING AROUSALS
i = i+1;
Ax(i) = plot_fig2_phaseangle(Fig, T_this, 'aw', 'phase_d0', pcfg);
Ax(i).Position = [0.5 0.22 0.18 0.3]+pcfg.margin;

% -------------------------------------------------------------------------
% Phase angle of CONT. SLEEP AROUSALS
i = i+1;
Ax(i) = plot_fig2_phaseangle(Fig, T_this, 'cs', 'phase_d0', pcfg);
Ax(i).Position = [0.67 0.22 0.18 0.3]+pcfg.margin;

% -------------------------------------------------------------------------
% Instantaneous amplitude of ISF
Ax(i) = plot_fig2_instamp(Fig, T_this, pcfg);
Ax(i).Position = [0.86 0.22 0.11 0.3]+pcfg.margin;

% -------------------------------------------------------------------------
% Panel labels
plot_fig2_panellabels(Fig, roi);

if strcmpi(roi, 'fz')
    Ax(2).Colormap = Ax(2).UserData.CMap;
    Ax(4).Colormap = Ax(4).UserData.CMap;
end
exportgraphics(Fig, sprintf('./figures/fig_prearousal_%s.png', roi), 'Resolution', 600)

end
