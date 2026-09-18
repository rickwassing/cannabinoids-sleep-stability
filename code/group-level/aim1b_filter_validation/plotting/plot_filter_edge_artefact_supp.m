function plot_filter_edge_artefact_supp(SIGMA, filtcfg, sig_test, sig_control, SIGA, SIGB, ANGA, ANGB, DANG, ford_mult, append_type)
% -------------------------------------------------------------------------
% Create supplementary figure showing the filter-edge-artefact validation
% method, the control/test signal traces, and the phase-difference
% distributions at increasing delays from the end of the bout.
% -------------------------------------------------------------------------
close all
Fig = figure();
Fig.Color = 'w';
Fig.Units = 'centimeters';
Fig.Position(3:4) = [18, 7];
clear Ax
ai = 0;
% -------------------------------------------------------------------------
% Plot the bars explaining the method
ai = ai+1;
Ax(ai) = axes('NextPlot', 'add');
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
XData = 1:2;
YData = [...
    length(sig_test.prepend), 4*filtcfg.order, 0, length(sig_test.append), 0, 0, 0, 0; ...
    0, 0, 0, 0, length(sig_control.prepend), 4*filtcfg.order, ford_mult*filtcfg.order, length(sig_control.append)];
bh = barh(XData, YData, 'stacked', 'BarWidth', 0.5);
plot(Ax(ai), [max(sum(YData'))-filtcfg.order, max(sum(YData'))], [1.2 1.2], '-k', 'LineWidth', 3)
text(Ax(ai), mean([max(sum(YData'))-filtcfg.order, max(sum(YData'))]), 1, 'filter length', 'HorizontalAlignment','center', 'VerticalAlignment','top', 'FontSize', 8)
plot(Ax(ai), [length(sig_test.prepend), 0], [0.75 0], ':k')
plot(Ax(ai), [length(sig_test.prepend)+4*filtcfg.order, max(sum(YData'))*0.74], [0.75 0], ':k')
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
bh(1).FaceColor = standard_colors('white');
bh(2).FaceColor = standard_colors('brick');
bh(3).FaceColor = standard_colors('brick').^0.6;
bh(4).FaceColor = standard_colors('white');
bh(5).FaceColor = standard_colors('white');
bh(6).FaceColor = standard_colors('cyan');
bh(7).FaceColor = standard_colors('cyan').^0.6;
bh(8).FaceColor = standard_colors('white');
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Axis props
Ax(ai).Box = 'off';
Ax(ai).Clipping = 'off';
Ax(ai).OuterPosition = [0/12 1-4/12 12/12 4/12] + [0 0 0.055 0];
Ax(ai).YLim = [0 2];
Ax(ai).YTick = 1:2;
Ax(ai).YTickLabel = {'test signal', 'control signal'};
Ax(ai).TickLength = [0, 0];
Ax(ai).XTick = [0.5*ford_mult*filtcfg.order, ford_mult*filtcfg.order, (4+ford_mult)*filtcfg.order, (4+2*ford_mult)*filtcfg.order, (4+2.5*ford_mult)*filtcfg.order];
Ax(ai).XTickLabel = {'prepend', sprintf('-%i', 4*filtcfg.order/SIGMA(1).srate), '0', sprintf('+%i', ford_mult*filtcfg.order/SIGMA(1).srate), 'append'};
Ax(ai).XLabel.String = 'time (s)';
Ax(ai).XAxisLocation = 'top';
Ax(ai).XAxis.Color = 'w';
Ax(ai).XAxis.Label.Color = 'k';
Ax(ai).XAxis.TickLabelColor = 'k';
Ax(ai).YAxis.Color = 'w';
Ax(ai).YAxis.TickLabelColor = 'k';
% -------------------------------------------------------------------------
% Plot the control and test signals
ai = ai+1;
Ax(ai) = axes('NextPlot', 'add', 'Layer', 'top');
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
idx_sel = 400;
patch(Ax(ai), 'XData', repmat(size(SIGA, 2), 1, 4)+[-500 0 0 -500], 'YData', [-8 -8 12 12], 'LineStyle', 'none', 'FaceColor', [0.92 0.93 0.95])
plot(1:size(SIGA, 2), SIGB(idx_sel, :).*4+5, '-', 'LineWidth', 1, 'Color', standard_colors('brick'))
plot(1:size(SIGA, 2), SIGA(idx_sel, :).*4+5, '-', 'LineWidth', 1, 'Color', standard_colors('cyan'))
plot(1:size(ANGB, 2), ANGB(idx_sel, :)-pi, '-', 'LineWidth', 1, 'Color', standard_colors('brick'))
plot(1:size(ANGA, 2), ANGA(idx_sel, :)-pi, '-', 'LineWidth', 1, 'Color', standard_colors('cyan'))
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Axis props
Ax(ai).Box = 'on';
Ax(ai).OuterPosition = [0/12 1-7/12 8/12 4/12]+[0.0425 0 0.045 -0.025];
Ax(ai).TickLength = [0, 0];
Ax(ai).XLim = [0, size(SIGA, 2)];
Ax(ai).XTick = [0, size(SIGA, 2)-500, size(SIGA, 2)];
Ax(ai).XTickLabel = {sprintf('-%i', 4*filtcfg.order/SIGMA(1).srate), '-50', '0'};
Ax(ai).YLim = [-8 12];
Ax(ai).YTick = [-pi 5];
Ax(ai).YTickLabel = {'angle', 'signal'};
% -------------------------------------------------------------------------
% Plot a Zoom-in of the control and test signals
ai = ai+1;
Ax(ai) = axes('NextPlot', 'add', 'Color', [0.92 0.93 0.95]);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
plot(1:size(SIGA, 2), SIGB(idx_sel, :).*4+5, '-', 'LineWidth', 1, 'Color', standard_colors('brick'))
plot(1:size(SIGA, 2), SIGA(idx_sel, :).*4+5, '-', 'LineWidth', 1, 'Color', standard_colors('cyan'))
plot(1:size(ANGB, 2), ANGB(idx_sel, :)-pi, '-', 'LineWidth', 1, 'Color', standard_colors('brick'))
plot(1:size(ANGA, 2), ANGA(idx_sel, :)-pi, '-', 'LineWidth', 1, 'Color', standard_colors('cyan'))
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Axis props
Ax(ai).Box = 'on';
Ax(ai).OuterPosition = [8/12 1-7/12 4/12 4/12]+[0.025 0 -0.04 -0.025];
Ax(ai).TickLength = [0, 0];
Ax(ai).XLim = [size(SIGA(idx_sel, :), 2)-500, size(SIGA(idx_sel, :), 2)];
Ax(ai).XTick = size(SIGA(idx_sel, :), 2)-500:50:size(SIGA(idx_sel, :), 2);
Ax(ai).XTickLabel = {'-50', '', '-40', '', '-30', '', '-20', '', '-10', '', '0'};
Ax(ai).XTickLabelRotation = 0;
Ax(ai).XGrid = 'on';
Ax(ai).YLim = Ax(ai-1).YLim;
Ax(ai).YTick = [];
% -------------------------------------------------------------------------
% Plot the difference in angle
ai = ai+1;
Ax(ai) = axes('NextPlot', 'add', 'Color', [0.92 0.93 0.95]);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
XData = ANGA(:, end-50:end)';
YData = DANG(:, end-50:end)';
XData = XData(:);
YData = YData(:);
[XData, idx_sort] = sort(XData);
YData = YData(idx_sort);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
scatter(Ax(ai), XData, YData, ...
    'SizeData', 1, ...
    'CData', [0 0 0], ...
    'Marker', 'o', ...
    'MarkerFaceAlpha', 0.01, ...
    'MarkerEdgeAlpha', 0.01);
EData = [];
EData.x = [];
EData.y = [];
EData.e = [];
for i = -pi:2*pi/18:pi-0.01
    idx_i = XData >=i & XData < i+2*pi/18;
    EData.x = [EData.x, i+2*pi/36];
    EData.y = [EData.y, circ_mean(YData(idx_i))];
    EData.e = [EData.e, circ_std(YData(idx_i))];
end
linepatch(Ax(ai), [-pi, pi], [0 0], 'EdgeColor', 'w', 'EdgeAlpha', 0.3)
errorbar(Ax(ai), EData.x, EData.y, EData.e, '.w', 'LineStyle', 'none', 'CapSize', 1)
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Axis props
Ax(ai).Box = 'on';
Ax(ai).OuterPosition = [8/12 1-12/12 4/12 5/12]+[0.055 0 -0.07 0];
Ax(ai).TickLength = [0, 0];
Ax(ai).XLim = [-pi, pi]+[-0.1, 0.1];
Ax(ai).YLim = [-pi, pi]+[-0.1, 0.1];
Ax(ai).XTick = [-pi, 0, pi];
Ax(ai).XTickLabel = {'-\pi', '0', '\pi'};
Ax(ai).YTick = [-pi, 0, pi];
Ax(ai).YTickLabel = {'-\pi', '0', '\pi'};
Ax(ai).XLabel.String = 'Phase angle (control)';
Ax(ai).YLabel.String = '\Deltaangle';
% -------------------------------------------------------------------------
% Legend for polar axis plots
ai = ai + 1;
Ax(ai) = polaraxes('NextPlot', 'add', 'Color', 'w');
Ax(ai).ThetaTick = [0, 90, 180, 270];
Ax(ai).ThetaTickLabel = {'0', '0.5\pi', '\pi', '1.5\pi'};
Ax(ai).FontSize = 8;
Ax(ai).RTick = [];
Ax(ai).OuterPosition = [(0)*1.175/12 1-12/12 1.2/12 5/12] + [0.01 0 0 0];
% -------------------------------------------------------------------------
% Plot the radial distributions
delay = [300 250 200 150 100 50 0];
for i = 1:length(delay)-1
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ai = ai + 1;
    Ax(ai) = polaraxes('NextPlot', 'add', 'Color', [0.92 0.93 0.95]);
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Ax(ai).OuterPosition = [(i)*1.175/12 1-12/12 1.2/12 5/12] + [0.01 0 0 0];
    Ax(ai).ThetaTick = [0, 90, 180, 270];
    Ax(ai).ThetaTickLabel = {};
    Ax(ai).RTick = [];
    Ax(ai).Title.String = sprintf('-%i to -%i s', delay(i)/10, delay(i+1)/10);
    Ax(ai).Title.FontSize = 8;
    Ax(ai).Title.FontWeight = 'normal';
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Theta = DANG(:, (end-delay(i)+1):(end-delay(i+1)));
    Theta = Theta(:);
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Mu = circ_mean(Theta);
    Sd = circ_std(Theta);
    fprintf('Mean (SD) angle difference at delay -%i to -%i is %.2f (%.2f).\n', delay(i), delay(i+1), Mu, Sd)
    h = polarhistogram(Ax(ai), Theta, 36, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', [0.5 0.5 0.5], 'Normalization', 'probability');
    polarplot(Ax(ai), [circ_mean(Theta), circ_mean(Theta)], [0 max(h.Values)], '-k', 'LineWidth', 1.5)
    Ax(ai).Subtitle.String = sprintf('%.2f° (%.2f°)', Mu, Sd);
    Ax(ai).Subtitle.FontSize = 8;
    Ax(ai).Subtitle.Units = 'normalized';
    Ax(ai).Subtitle.Position(2) = -0.30;
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Ax(ai).RLim = [0, max(h.Values)];
end
% -------------------------------------------------------------------------
% Save image
exportgraphics(Fig, sprintf('figures/supp_2_append-%s.png', append_type), 'Resolution', 300)
disp('Done saving')

end
