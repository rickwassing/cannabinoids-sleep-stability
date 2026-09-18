function [Ax, BData] = plot_fig3_mdlhists(desmat, fld, BinEdges, orient)

Ax = axes();

BinCenters = BinEdges(1:end-1)+mean(diff(BinEdges))/2;
clear BData
BData.pbo = histcounts(desmat.(fld)(strcmpi(desmat.cond, 'placebo')), 'BinEdges', BinEdges)';
BData.etc = histcounts(desmat.(fld)(strcmpi(desmat.cond, 'etc120')), 'BinEdges', BinEdges)';
switch orient
    case 'h'
        h = barh(Ax, BinCenters, [BData.pbo, BData.etc], ...
            'BarWidth', 0.7, ...
            'GroupWidth', 1);
    case 'v'
        h = bar(Ax, BinCenters, [BData.pbo, BData.etc], ...
            'BarWidth', 0.7, ...
            'GroupWidth', 1);
end

h(1).EdgeColor = css_standard_colors('pbo');
h(1).FaceColor = css_standard_colors('pbo');
h(2).EdgeColor = css_standard_colors('etc');
h(2).FaceColor = css_standard_colors('etc');

Ax.Color = [0.96 0.97 0.99];
Ax.Box = 'on';
Ax.FontSize = 8;
Ax.TickLength = [0, 0];
switch orient
    case 'h'
        Ax.XTick = Ax.XTick(end);
        Ax.YTick = 0:60:420;
        Ax.YTickLabel = [];
        Ax.YLim = [BinEdges(1), BinEdges(end)];
        Ax.YGrid = 'on';
    case 'v'
        Ax.YTick = [];
        Ax.XTick = BinEdges;
        Ax.XTickLabel = [];
        Ax.XLim = [BinEdges(1), BinEdges(end)];
end
end