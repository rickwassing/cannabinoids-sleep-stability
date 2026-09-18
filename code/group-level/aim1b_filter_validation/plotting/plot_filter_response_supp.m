function plot_filter_response_supp(SIGMA, filtcfg, append_type)
% -------------------------------------------------------------------------
% Plot the filter magnitude/phase response and export as a supplementary
% figure.
% -------------------------------------------------------------------------
close all
pop_firws(SIGMA(1), ...
    'fcutoff', filtcfg.cutoff+filtcfg.adj, ...
    'ftype', 'bandpass', ...
    'wtype', filtcfg.wintype, ...
    'warg', filtcfg.warg, ...
    'forder', filtcfg.order, ...
    'plotfresp', true, ...
    'minphase', 0);
Parent = gcf;
PhaseAx = Parent.Children(1);
MagAx = Parent.Children(2);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Fig = figure();
Fig.Units = 'centimeters';
Fig.Position(3:4) = [17 8];
Ax(1) = copyobj(MagAx, Fig);
Ax(1).OuterPosition = [0 0 0.5 1];
Ax(1).XLim = [0 0.1];
Ax(1).YLim = [-70 1];
Ax(1).YTick = -60:10:10;
Ax(1).XTick = sort([0.1, ...
    MagAx.Children.XData(find(MagAx.Children.YData > -20, 1, 'first')), ...
    MagAx.Children.XData(find(MagAx.Children.YData == max(MagAx.Children.YData), 1, 'first')), ...
    MagAx.Children.XData(find(MagAx.Children.YData > -20, 1, 'last')), ...
    ]);
Ax(1).XTickLabelRotation = 90;
Ax(1).FontSize = 10;
Ax(2) = copyobj(PhaseAx, Fig);
Ax(2).OuterPosition = [0.5 0 0.5 1];
Ax(2).XLim = [0 0.1];
Ax(2).YLim = [-pi pi];
Ax(2).YTick = -pi:0.5*pi:pi;
Ax(2).XTick = sort([0.1, ...
    MagAx.Children.XData(find(MagAx.Children.YData > -20, 1, 'first')), ...
    MagAx.Children.XData(find(MagAx.Children.YData == max(MagAx.Children.YData), 1, 'first')), ...
    MagAx.Children.XData(find(MagAx.Children.YData > -20, 1, 'last')), ...
    ]);
Ax(2).XTickLabelRotation = 90;
Ax(2).YLim = [-0.5 3.6];
Ax(2).YTickLabel = {'-\pi', '', '0', '', '\pi'};
Ax(2).FontSize = 10;
Ax(1).Children(1).Marker = '.';
Ax(2).Children(1).Marker = '.';
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fname = struct();
fname.fcut = strrep(sprintf('%.4f%.4f',filtcfg.cutoff), '0.', 'x');
fname.tbw = sprintf('1x%i', 1/filtcfg.transbw);
fname.rdev = strrep(sprintf('%.4f',filtcfg.rippledev), '0.', 'x');
exportgraphics(Fig, sprintf('./figures/filtresponse_ftype-%s_forder-%i_fcut-%s_tbw-%s_rdev-%s_atype-%s.png', filtcfg.wintype, filtcfg.order, fname.fcut, fname.tbw, fname.rdev, strrep(append_type, '_', '')), 'Resolution', 300)

end
