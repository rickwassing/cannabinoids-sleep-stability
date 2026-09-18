function plot_supp_phase_coupling(T_this, pcfg)
% -------------------------------------------------------------------------
% Supplementary figure: phase-angle polar histograms of awakening /
% continued-sleep arousals at increasing delays before the sigma ISF
% zero-crossing (phase_d0, phase_d4, phase_d8, phase_d12).
% -------------------------------------------------------------------------
Fig = figure('Color', 'w');
Fig.Units = 'centimeters';
Fig.Position = [1 12 8.5 16].*1;

clear Ax
i = 0;

i = i+1;
Ax(i) = plot_fig2_phaseangle(Fig, T_this, 'aw', 'phase_d0', pcfg);
Ax(i).OuterPosition = [0 0.75 0.5 0.25];

i = i+1;
Ax(i) = plot_fig2_phaseangle(Fig, T_this, 'cs', 'phase_d0', pcfg);
Ax(i).OuterPosition = [0.5 0.75 0.5 0.25];

i = i+1;
Ax(i) = plot_fig2_phaseangle(Fig, T_this, 'aw', 'phase_d4', pcfg);
Ax(i).OuterPosition = [0 0.5 0.5 0.25];

i = i+1;
Ax(i) = plot_fig2_phaseangle(Fig, T_this, 'cs', 'phase_d4', pcfg);
Ax(i).OuterPosition = [0.5 0.5 0.5 0.25];

i = i+1;
Ax(i) = plot_fig2_phaseangle(Fig, T_this, 'aw', 'phase_d8', pcfg);
Ax(i).OuterPosition = [0 0.25 0.5 0.25];

i = i+1;
Ax(i) = plot_fig2_phaseangle(Fig, T_this, 'cs', 'phase_d8', pcfg);
Ax(i).OuterPosition = [0.5 0.25 0.5 0.25];

i = i+1;
Ax(i) = plot_fig2_phaseangle(Fig, T_this, 'aw', 'phase_d12', pcfg);
Ax(i).OuterPosition = [0 0 0.5 0.25];

i = i+1;
Ax(i) = plot_fig2_phaseangle(Fig, T_this, 'cs', 'phase_d12', pcfg);
Ax(i).OuterPosition = [0.5 0 0.5 0.25];

exportgraphics(Fig, './figures/supp_2b_phase-coupling_fz.png', 'Resolution', 600)

end
