
Fig = figure('Position', [10, 300, 560, 420]);

Ax = axes(Fig);
Ax.Layer = 'top'; 
Ax.Box = 'on';
Ax.TickDir = 'out';
Ax.NextPlot = 'add';

XData = -90:0.5:-0.5;
YData = 1:30;
AvCData = squeeze(mean(MdlData.logpow.awake.tstat, 1))';

SigData = squeeze(sum(MdlData.logpow.awake.pval < 0.05, 1))';

imagesc('XData', XData, 'YData', YData, 'CData', AvCData)
contour(XData, YData, SigData, 1, 'k');

Ax.XLim = [-90.5, 0];
Ax.YLim = [0.5; 30.5];
Ax.CLim = [-3.1, 3.1];


Fig = figure('Position', [450, 300, 560, 420]);

Ax = axes(Fig);
Ax.Layer = 'top'; 
Ax.Box = 'on';
Ax.TickDir = 'out';
Ax.NextPlot = 'add';

XData = -90:0.5:-0.5;
YData = 1:30;
AvCData = squeeze(mean(MdlData.logpow.cond.tstat, 1))';

SigData = squeeze(sum(MdlData.logpow.cond.pval < 0.05, 1))';

imagesc('XData', XData, 'YData', YData, 'CData', AvCData)
contour(XData, YData, SigData, 1, 'k');

Ax.XLim = [-90.5, 0];
Ax.YLim = [0.5; 30.5];
Ax.CLim = [-3.1, 3.1];


Fig = figure('Position', [900, 300, 560, 420]);

Ax = axes(Fig);
Ax.Layer = 'top'; 
Ax.Box = 'on';
Ax.TickDir = 'out';
Ax.NextPlot = 'add';

XData = -90:0.5:-0.5;
YData = 1:30;
AvCData = squeeze(mean(MdlData.logpow.intx.tstat, 1))';

SigData = squeeze(sum(MdlData.logpow.intx.pval < 0.05, 1))';

imagesc('XData', XData, 'YData', YData, 'CData', AvCData)
contour(XData, YData, SigData, 1, 'k');

Ax.XLim = [-90.5, 0];
Ax.YLim = [0.5; 30.5];
Ax.CLim = [-3.1, 3.1];

%%

t = -11.4:-4;
f = 14:15;

times = -90:0.5:30;
specfreqs = 1:1:30;

idx_t = times >= min(t) & times <= max(t);
idx_f = specfreqs >= min(f) & specfreqs <= max(f);

idx_e = strcmpi(Cond.cs, 'etc120');
Y_cs_plc = log10(squeeze(mean(mean(mean(CData.cs(:, idx_t, idx_f, ~idx_e), 3, 'omitnan'), 2, 'omitnan'), 1, 'omitnan')));
Y_cs_etc = log10(squeeze(mean(mean(mean(CData.cs(:, idx_t, idx_f, idx_e), 3, 'omitnan'), 2, 'omitnan'), 1, 'omitnan')));

idx_e = strcmpi(Cond.aw, 'etc120');
Y_aw_plc = log10(squeeze(mean(mean(mean(CData.aw(:, idx_t, idx_f, ~idx_e), 3, 'omitnan'), 2, 'omitnan'), 1, 'omitnan')));
Y_aw_etc = log10(squeeze(mean(mean(mean(CData.aw(:, idx_t, idx_f, idx_e), 3, 'omitnan'), 2, 'omitnan'), 1, 'omitnan')));

G = [...
    repmat({'cs_plc'}, 1, length(Y_cs_plc)), ...
    repmat({'cs_etc'}, 1, length(Y_cs_etc)), ...
    repmat({'aw_plc'}, 1, length(Y_aw_plc)), ...
    repmat({'aw_etc'}, 1, length(Y_aw_etc)), ...
    ];

Fig = figure;
boxplot([Y_cs_plc; Y_cs_etc; Y_aw_plc; Y_aw_etc], G)

