function Ax = plot_fig3_ampdur(desmat, fld)

Ax = axes('NextPlot', 'add');
Ax.Color = [0.96 0.97 0.99];
Ax.Box = 'on';
Ax.YGrid = 'on';
Ax.FontSize = 8;
Ax.TickLength = [0, 0];
Ax.XLim = [0.33, 4.67];
Ax.XTick = [1.5, 3.5];
Ax.XTickLabel = {'NREM', 'TR'};
switch fld
    case 'amplitude'
        Ax.YLabel.String = 'Amp. (dB)';
        Ax.YLim = [0.67 1.1];
        Ax.YTick = [0.84, 0.95];
    case 'duration'
        Ax.YLabel.String = 'Dur. (s)';
        Ax.YLim = [20 35];
        Ax.YTick = [25, 30];
        Ax.YTickLabel = Ax.YTick;
end

Ax.YLabel.FontSize = 10;

clear YData

YData.nrem.pbo.d = desmat.(fld)(strcmpi(desmat.cond, 'placebo') & strcmpi(desmat.stage, 'nrem'));
YData.rem.pbo.d = desmat.(fld)(strcmpi(desmat.cond, 'placebo') & strcmpi(desmat.stage, 'rem'));
YData.nrem.etc.d = desmat.(fld)(strcmpi(desmat.cond, 'etc120') & strcmpi(desmat.stage, 'nrem'));
YData.rem.etc.d = desmat.(fld)(strcmpi(desmat.cond, 'etc120') & strcmpi(desmat.stage, 'rem'));
YData.nrem.pbo.mu = mean(YData.nrem.pbo.d);
YData.rem.pbo.mu = mean(YData.rem.pbo.d);
YData.nrem.etc.mu = mean(YData.nrem.etc.d);
YData.rem.etc.mu = mean(YData.rem.etc.d);

YData.nrem.pbo.e = tinv(0.975, length(YData.nrem.pbo.d)-1).*(std(YData.nrem.pbo.d)./sqrt(length(YData.nrem.pbo.d)));
YData.rem.pbo.e = tinv(0.975, length(YData.rem.pbo.d)-1).*(std(YData.rem.pbo.d)./sqrt(length(YData.rem.pbo.d)));
YData.nrem.etc.e = tinv(0.975, length(YData.nrem.etc.d)-1).*(std(YData.nrem.etc.d)./sqrt(length(YData.nrem.etc.d)));
YData.rem.etc.e = tinv(0.975, length(YData.rem.etc.d)-1).*(std(YData.rem.etc.d)./sqrt(length(YData.rem.etc.d)));

errorbar(Ax, [1 3], [YData.nrem.pbo.mu, YData.rem.pbo.mu], [YData.nrem.pbo.e, YData.rem.pbo.e], 'o', ...
    'LineStyle', 'none', ...
    'Color', 'k', ...
    'MarkerSize', 3, ...
    'MarkerFaceColor', css_standard_colors('pbo'), ...
    'MarkerEdgeColor', css_standard_colors('pbo'));

errorbar(Ax, [2 4], [YData.nrem.etc.mu, YData.rem.etc.mu], [YData.nrem.etc.e, YData.rem.etc.e], 'o', ...
    'LineStyle', 'none', ...
    'Color', 'k', ...
    'MarkerSize', 3, ...
    'MarkerFaceColor', css_standard_colors('etc'), ...
    'MarkerEdgeColor', css_standard_colors('etc'));

plot(Ax, [0.6, 2.33], [Ax.YLim(1), Ax.YLim(1)], '-', ...
    'Color', css_standard_colors('blue'), ...
    'LineWidth', 3)

plot(Ax, [2.6, 4.33], [Ax.YLim(1), Ax.YLim(1)], '-', ...
    'Color', css_standard_colors('red'), ...
    'LineWidth', 3)

end