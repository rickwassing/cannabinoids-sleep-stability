function [Ax] = plot_fig2_eegtrace(Fig, EEG, fld, pcfg)

% Get indexes
idx_ar = round(EEG.event(pcfg.idx.aro.(fld)).latency);
idx_x = idx_ar+pcfg.xlim(1)*EEG.srate:idx_ar+pcfg.xlim(2)*EEG.srate;

% Create axes
Ax = axes(Fig, ...
    'NextPlot', 'add');

% Axis props
Ax.XAxis.Visible = 'off';
Ax.YAxis.Color = 'w';
Ax.YLim = pcfg.ylim_eeg;
Ax.XLim = [idx_x(1), idx_x(end)]./EEG.srate;
Ax.YTick = [];
Ax.YLabel.String = 'EEG';
Ax.YLabel.FontSize = 10;
Ax.YLabel.Color = 'k';

% Plot EEG data
XData = EEG.times(idx_x);
YData = EEG.data(1, idx_x);
plot(Ax, XData, YData, '-k', ...
    'LineWidth', 0.25)

% Plot arousal marker
XData = EEG.times(idx_ar);
plot(Ax, [XData, XData], pcfg.ylim_eeg, ':k')

% Write text
switch fld
    case 'aw'
        str = ' N2 → Arousal → Wake ';
        clr = 'aw';
    case 'cs'
        str = ' N2 → Arousal → N2 ';
        clr = 'cs';
end
text(Ax, XData-7.5, pcfg.ylim_eeg(2), str, ...
    'FontSize', 8, ...
    'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'top', ...
    'Color', css_standard_colors(clr))

end