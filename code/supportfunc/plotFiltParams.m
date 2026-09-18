function [Fig, fname] = plotFiltParams(SIGMA, filtcfg, const)

% Sampling rate
fs = SIGMA(1).srate;

% Normalize cutoffs
wn = filtcfg.cutoff / (fs/2);

% Design FIR bandpass filter
b = fir1(filtcfg.order, wn, 'bandpass', kaiser(filtcfg.order+1, filtcfg.warg));

% Plot filter coefficients
close all
pop_firws(SIGMA(1), ...
    'fcutoff', filtcfg.cutoff, ...
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
Fig.Position(3:4) = [17 6];

clear Ax

Ax(1) = copyobj(MagAx, Fig);
Ax(1).OuterPosition = [0 0 0.33 1];
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
Ax(2).OuterPosition = [0.33 0 0.33 1];
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

Ax(3) = axes(Fig);
Ax(3).OuterPosition = [0.67 0 0.33 1];
Ax(3).Position(2) = Ax(2).Position(2);
Ax(3).Position(4) = Ax(2).Position(4);
plot(0:filtcfg.order, b, 'Color', [0.8500 0.3250 0.0980]);
Ax(3).XTick = [0 130 315 500 630];
Ax(3).YLim = [-0.002 0.008];
Ax(3).FontSize = 10;
grid on
xlabel('Coefficient index (n)');
ylabel('Coefficients');
Ax(3).Title.String ='Filter Coefficients';
Ax(3).Title.Position(1) = 385;

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fname = struct();
fname.fcut = strrep(sprintf('%.4f%.4f',filtcfg.cutoff), '0.', 'x');
fname.tbw = sprintf('1x%i', 1/filtcfg.transbw);
fname.rdev = strrep(sprintf('%.4f',filtcfg.rippledev), '0.', 'x');
fname.warg = sprintf('%i', filtcfg.warg);
fname = sprintf('./figures/filtresponse_ftype-%s_forder-%i_fcut-%s_tbw-%s_rdev-%s_warg-%s_atype-%s.png', filtcfg.wintype, filtcfg.order, fname.fcut, fname.tbw, fname.rdev, fname.warg, strrep(const.append_type, '_', ''));

end