function plot_fig3_full(SIG, HYP, bouts, chans, H, desmat, chans_sel)
% -------------------------------------------------------------------------
% Assemble Figure 3: hypnogram with selected bouts, average sigma trace at
% REM transition, amplitude/duration bar plots, channel topoplot, and the
% linear-model panels with their marginal histograms.
% -------------------------------------------------------------------------
clear pcfg
pcfg.markersize = 5;

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
        'FaceColor', standard_colors('blue'), ...
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

if strcmpi(chans_sel, 'fz')
    exportgraphics(Fig, './figures/figure3_fz.png', 'Resolution', 1200)
else
    exportgraphics(Fig, './figures/figure3_pz.png', 'Resolution', 1200)
end

end
