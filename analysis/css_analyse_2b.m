%function css_analyse_2
% -------------------------------------------------------------------------
% Use the the 130-second pre-arousal bouts of continuous N2 sleep to
% determine the phase angle, amplitude (hilbert?) between state-shift
% arousals and continued sleep arousals; and whether CBN modulates this
% - State-shift arousals occur after the sigma ISF peak, and
%   continued-sleep arousals occur prior to the peak.
% - State-shift arousals occur at higher ISF amplitudes (hilbert?)
% - CBN increases the number of arousals to after the peak, or at higher
%   ISF amplitudes
% -------------------------------------------------------------------------
% Define the set of parietal channels to focus on
clc
clear
close all
chanlocs = template_to_chanlocs(which('GSN-HydroCel-257.sfp'));
[~, system] = ishdeeg({chanlocs.labels});
incl = hdeeg_scalpchannels(system);
chanlocs = chanlocs(ismember({chanlocs.labels}, incl));
chanlocs = channel_clusters(chanlocs, 'mff');
idx_chan = [77, 78, 79, 85, 86, 87, 92, 93, 94, 99, 104, 110, 111, 112, 113, 120, 121, 122];
chanlabels = {chanlocs(idx_chan).labels};
append_type = 'zeros';
ford_mult = 2;
% -------------------------------------------------------------------------
% Create supplementary figure showing which channels were selected
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Load ISF data
type = 'norm';
Files = dir(sprintf('derivatives/EEG-output-fstlvl/sub-*/ses-*/sub-*%ssigma*interp_fstlvl.mat', type));
ISF = [];
for i = 1:length(Files)
    if i == 1
        ISF = LoadDataset(fullfile(Files(i).folder, Files(i).name), 'matrix');
    else
        ISF(i) = LoadDataset(fullfile(Files(i).folder, Files(i).name), 'matrix'); %#ok<SAGROW>
    end
end
YData = arrayfun(@(s) s.features(2).data, ISF, 'UniformOutput', false);
YData = cat(2, YData{:});
YData = mean(YData, 2);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
Fig = figure();
Fig.Units = 'centimeters';
Fig.Position(3:4) = [8 8];
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Plot as topoplot
load('colormap_batlow.mat')
topoplot(YData, chanlocs, ...
    'colormap', batlow, ...
    'hlinewidth', 1, ...
    'maplimits', [0, 2], ...
    'numcontour', 0, ...
    'conv', 'on', ...
    'whitebk', 'on', ...
    'emarker2', {idx_chan, '.', 'w', 15, 1});
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
exportgraphics(Fig, 'figures/supp_selparietalchans.png', 'Resolution', 600)
% -------------------------------------------------------------------------
% Load segmented pre-arousal bouts
Files = dir('derivatives/EEG-segmented/sub-*/ses-*/sub-*_desc-sigmanremarobout*.set');
clear SIGMA
for i = 1:length(Files)
    SIGMA(i) = LoadDataset(fullfile(Files(i).folder, Files(i).name), 'all'); %#ok<SAGROW>
    SIGMA(i) = pop_select(SIGMA(i), 'channel', idx_chan); %#ok<SAGROW>
end
% -------------------------------------------------------------------------
% Emperical values for the filter cutoffs
Files = dir('derivatives/EEG-output-fstlvl/sub-*/ses-*/sub-*_desc-a1cabssigma_fstlvl.mat');
filtcfg = struct();
filtcfg.cutoff = [];
for i = 1:length(Files)
    tmp = LoadDataset(fullfile(Files(i).folder, Files(i).name), 'matrix');
    bw = [...
        tmp.features(1).data(idx_chan) - tmp.features(3).data(idx_chan)/2, ...
        tmp.features(1).data(idx_chan) + tmp.features(3).data(idx_chan)/2];
    filtcfg.cutoff = [filtcfg.cutoff; mean(bw)];
end
% Get filter cutoffs and filter order
% filtcfg.cutoff = mean(filtcfg.cutoff, 1); % 0.0138 0.0259
% filtcfg.wintype = 'kaiser';
% filtcfg.transbw = 1/50;
% filtcfg.rippledev = 0.005;
% filtcfg.adj = [0 0];%[0.005 -0.001];
% filtcfg.order = pop_firwsord(filtcfg.wintype, SIGMA(1).srate, filtcfg.transbw, filtcfg.rippledev); % transition bandwidth of 0.01 Hz
% filtcfg.order

filtcfg.cutoff = mean(filtcfg.cutoff, 1); % 0.0138 0.0259
filtcfg.wintype = 'kaiser';
filtcfg.transbw = 1/50;
filtcfg.rippledev = 0.05;
filtcfg.warg = 2;
filtcfg.adj = [0 0];%[0.005 -0.001];
filtcfg.order = pop_firwsord(filtcfg.wintype, SIGMA(1).srate, filtcfg.transbw, filtcfg.rippledev); % transition bandwidth of 0.01 Hz
% filtcfg.order = 1000;
filtcfg.order

% -------------------------------------------------------------------------
% Plot the filter response
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
Fig.Position(3:4) = [24 12];
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
Ax(1).FontSize = 16;
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
Ax(2).FontSize = 16;
Ax(1).Children(1).Marker = '.';
Ax(2).Children(1).Marker = '.';
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fname = struct();
fname.fcut = strrep(sprintf('%.4f%.4f',filtcfg.cutoff), '0.', 'x');
fname.tbw = sprintf('1x%i', 1/filtcfg.transbw);
fname.rdev = strrep(sprintf('%.4f',filtcfg.rippledev), '0.', 'x');
exportgraphics(Fig, sprintf('./figures/filtresponse_ftype-%s_forder-%i_fcut-%s_tbw-%s_rdev-%s_atype-%s.png', filtcfg.wintype, filtcfg.order, fname.fcut, fname.tbw, fname.rdev, strrep(append_type, '_', '')), 'Resolution', 300)

% -------------------------------------------------------------------------
% Crop the signal by 31.5 seconds from the end-boundary (30 seconds were
% added to include the arousal, and the wavelett at 11 Hz had a duration of
% 1.45 seconds (16 cycles * 1/11 seconds per cycle).
crop_aro = [0, -31.5].*SIGMA(1).srate;
% -------------------------------------------------------------------------
% Init output table
T = [];
Pr_aro = struct();
Pr_aro.cs = [];
Pr_aro.wk = [];
% -------------------------------------------------------------------------
% init checking vars
times = [];
inspect = struct();
inspect.angle = struct();
inspect.raw = struct();
inspect.pred = struct();
inspect.filt = struct();
inspect.angle.x = []; inspect.angle.y = []; inspect.angle.m = [];
inspect.raw.x = []; inspect.raw.y = [];
inspect.pred.x = []; inspect.pred.y = [];
inspect.filt.x = []; inspect.filt.y = [];
% -------------------------------------------------------------------------
% Loop over files and extract bouts
cbouts = struct();
cbouts.blen = [];
cbouts.total = 0;
cbouts.exec = 0;
cbouts.tooshort = 0;
cbouts.noaro = 0;
cbouts.precaro = 0;
for i = 1:length(SIGMA)
    % Get indices of the selected arousal bouts
    bouts = find(strcmpi({SIGMA(i).event.type}, 'boundary'));
    bouts = ceil([1, [SIGMA(i).event(bouts).latency], SIGMA(i).pnts]);
    % Loop over bouts
    for b = 1:length(bouts)-1
        cbouts.total = cbouts.total+1;
        % ... skip if the bout is shorter than the filter order
        switch append_type
            case 'none'
                boutlength = (bouts(b+1)-bouts(b))+crop_aro(2);
            otherwise
                boutlength = (bouts(b+1)-bouts(b))+2*ford_mult*filtcfg.order;
        end
        cbouts.blen = [cbouts.blen, boutlength];
        % Check if bout is long enough
        if boutlength < 2*filtcfg.order+1
            cbouts.tooshort = cbouts.tooshort+1;
            continue
        end
        % Get the index of the selected event
        idx_evt = find([SIGMA(i).event.latency] > bouts(b+1)-30 * SIGMA(i).srate & [SIGMA(i).event.latency] < bouts(b+1)-29.8 * SIGMA(i).srate);
        idx_evt = idx_evt(contains({SIGMA(i).event(idx_evt).type}, 'arousal'));
        if isempty(idx_evt)
            cbouts.noaro = cbouts.noaro+1;
            continue
        end
        if length(idx_evt) ~= 1
            error('found more than one selected arousal!')
        end
        % Create filtered signal
        sigaro = pop_select(SIGMA(i), 'point', [bouts(b), bouts(b+1)]);
        sig = pop_select(SIGMA(i), 'point', [bouts(b), bouts(b+1)]+crop_aro);
        % Check if an arousal occurred prior to the selected arousal
        if any(contains({sig.event.type}, 'arousal'))
            if any([sig.event(find(contains({sig.event.type}, 'arousal'))).duration] > 15.*sig.srate)
                cbouts.precaro = cbouts.precaro+1;
                continue
            end
        end
        cbouts.exec = cbouts.exec+1;
        % Average across selected channels and de-mean
        sigaro.data = detrend(mean(sigaro.data, 1, 'omitnan')', 0)';
        sigaro.nbchan = size(sig.data, 1);
        sigaro.times = sigaro.times./1000;
        sig.data = detrend(sig.data', 0)';
        sig.nbchan = size(sig.data, 1);
        sig.times = sig.times./1000;
        % Apply one of the following methods to minimize the edge artefact
        switch append_type
            case 'flip'
                % Append a flipped copy of itself at the end of the signal
                if sig.pnts < ford_mult*filtcfg.order
                    sig.prepend = zeros(sig.nbchan, ford_mult*filtcfg.order);
                    sig.append = zeros(sig.nbchan, ford_mult*filtcfg.order);
                    sig.prepend(:, (end-sig.pnts)+1:end) = doubleflip(sig.data', 'none')'; %-1.*flip(sig.data, 2);
                    sig.append(:, 1:sig.pnts) = doubleflip(sig.data', 'none')'; %-1.*flip(sig.data, 2);
                else
                    sig.prepend = doubleflip(sig.data(:, 1:ford_mult*filtcfg.order)', 'none');
                    sig.append = doubleflip(sig.data(:, end-ford_mult*filtcfg.order+1:end)', 'none');
                end
            case 'zeros'
                % Append with zeros
                sig.prepend = zeros(sig.nbchan, ford_mult*filtcfg.order);
                sig.append = zeros(sig.nbchan, ford_mult*filtcfg.order);
            case 'autoreg'
                % Autoregressive prediction
                if sig.pnts < ford_mult*filtcfg.order
                    sig.prepend = flip(ar_pred(flip(sig.data, 2)', ford_mult*filtcfg.order, sig.pnts)', 2);
                    sig.append = ar_pred(sig.data', ford_mult*filtcfg.order, sig.pnts)';
                else
                    sig.prepend = flip(ar_pred(flip(sig.data, 2)', ford_mult*filtcfg.order)', 2);
                    sig.append = ar_pred(sig.data', ford_mult*filtcfg.order)';
                end
            case 'none'
                sig.prepend = [];
                sig.append = [];
            otherwise
                error('Specify an append type.')
        end
        % Append data
        sig.data = [sig.prepend, sig.data, sig.append];
        sig.pnts = size(sig.data, 2);
        sig.xmin = 0;
        sig.xmax = sig.pnts/sig.srate;
        sig.times = linspace(sig.xmin, sig.xmax, sig.pnts);
        % Filter data within ISF freq range
        fsig = pop_firws(sig, ...
            'fcutoff', filtcfg.cutoff+filtcfg.adj, ...
            'ftype', 'bandpass', ...
            'wtype', filtcfg.wintype, ...
            'warg', filtcfg.warg, ...
            'forder', filtcfg.order, ...
            'plotfresp', false, ...
            'minphase', 0);
        % Calculate Hilbert on filtered signal
        HilbAngle = hilbert(detrend(fsig.data', 0));
        HilbAngle = angle(HilbAngle); % Get phase-angle
        MuHilbAngle = circ_mean(HilbAngle');
        idx_aro_onset = length(MuHilbAngle) - length(sig.append);
        % TODO unwarp
        inspect.angle.x = [inspect.angle.x, {fsig.times - fsig.times(idx_aro_onset)}]; inspect.angle.y = [inspect.angle.y, {MuHilbAngle}]; inspect.angle.m = [inspect.angle.m, MuHilbAngle(idx_aro_onset)];
        inspect.raw.x = [inspect.raw.x, {sigaro.times - fsig.times(idx_aro_onset)}]; inspect.raw.y = [inspect.raw.y, {sigaro.data}];
        inspect.pred.x = [inspect.pred.x, {sig.times - fsig.times(idx_aro_onset)}]; inspect.pred.y = [inspect.pred.y, {sig.data}];
        inspect.filt.x = [inspect.filt.x, {fsig.times - fsig.times(idx_aro_onset)}]; inspect.filt.y = [inspect.filt.y, {fsig.data}];
        % Create new row for output table
        kv = filename2struct(SIGMA(i).setname);
        tmp = table();
        tmp.sub = {kv.sub};
        tmp.ses = {kv.ses};
        tmp.bout = b;
        tmp.aro_type = {SIGMA(i).event(idx_evt).type};
        tmp.is_awakening = SIGMA(i).event(idx_evt).is_awakening;
        tmp.stage = SIGMA(i).event(idx_evt).stage;
        tmp.next_stage = SIGMA(i).event(idx_evt).next_stage;
        tmp.duration = SIGMA(i).event(idx_evt).duration;
        tmp.origlatency = SIGMA(i).event(idx_evt).origlatency;
        tmp.phase_d30 = MuHilbAngle(idx_aro_onset-29.*fsig.srate);
        tmp.phase_d25 = MuHilbAngle(idx_aro_onset-24.*fsig.srate);
        tmp.phase_d20 = MuHilbAngle(idx_aro_onset-19.*fsig.srate);
        tmp.phase_d15 = MuHilbAngle(idx_aro_onset-14.*fsig.srate);
        tmp.phase_d10 = MuHilbAngle(idx_aro_onset-9.*fsig.srate);
        tmp.phase_d5 = MuHilbAngle(idx_aro_onset-4.*fsig.srate);
        tmp.phase_d4 = MuHilbAngle(idx_aro_onset-3.*fsig.srate);
        tmp.phase_d3 = MuHilbAngle(idx_aro_onset-2.*fsig.srate);
        tmp.phase_d2 = MuHilbAngle(idx_aro_onset-1.*fsig.srate);
        tmp.phase_d1 = MuHilbAngle(idx_aro_onset);
        tmp.phase_cell = {HilbAngle(idx_aro_onset-0.*fsig.srate, :)};
        % store the row in the output table
        if isempty(T)
            T = tmp;
        else
            T = [T; tmp]; %#ok<AGROW>
        end
        % As a control, get box-car timeseries for other arousals in the bout
        x = events2timeseries(fsig.event(contains({fsig.event.type}, 'arousal')), fsig.xmin, fsig.xmax, fsig.srate);
        x = x(end-30*fsig.srate+1:end);
        if tmp.is_awakening
            Pr_aro.wk = [Pr_aro.wk; x];
        else
            Pr_aro.cs = [Pr_aro.cs; x];
        end
    end
end

%%

session = {'placebo'};
aro_type = {'arousalemg'};

idx_sel = find(ismember(T.aro_type, aro_type) & ismember(T.ses, session));

[~, idx_sort] = sort(inspect.angle.m(idx_sel));
idx_sel = idx_sel(idx_sort);

XTick = -200:50:200;

nplt = 15;
nfigs = ceil(length(idx_sel)/nplt);

for f = 1:nfigs
    close all
    Fig = figure('Position', [680 100 630 900]);
    Ax = axes(Fig, 'NextPlot', 'add');
    Ax.TickLength = [0, 0];
    Ax.XTick = XTick;
    Ax.XTickLabel = repmat({''}, 1, length(XTick));
    Ax.XTickLabel(XTick == 0) = {'0'};
    Ax.XGrid = 'on';
    Ax.Color = [0.95 0.96 0.98];
    Ax.YTick = [];
    Ax.YLim = [25 nplt*50+25];
    Ax.Position = [0 0 0.65 1];
    plot(Ax, [0 0], Ax.YLim, '-k')
    if max((1:nplt)+(f-1)*nplt) > length(idx_sel)
        this_set = idx_sel(1+(f-1)*nplt:end);
    else
        this_set = idx_sel((1:nplt)+(f-1)*nplt);
    end
    cnt = 0;
    for i = asrow(this_set)
        cnt = cnt+1;
    %     plot(Ax, inspect.raw.x{i}, inspect.raw.y{i}'+cnt*50, '-', 'Color', [0.67 0.67 0.67])
        linepatch(Ax, inspect.pred.x{i}, inspect.pred.y{i}'+cnt*50, 'LineWidth', 1, 'EdgeAlpha', 0.01)
        linepatch(Ax, inspect.filt.x{i}, inspect.filt.y{i}'*10+cnt*50+10, 'LineWidth', 1, 'EdgeAlpha', 0.25)
        plot(Ax, inspect.angle.x{i}, inspect.angle.y{i}'*2+cnt*50-4*pi, '-', 'LineWidth', 2, 'Color', standard_colors('red'))
        plot(Ax, 0, inspect.angle.m(i)*2+cnt*50-4*pi, 'o', 'MarkerFaceColor', standard_colors('red'), 'MarkerEdgeColor', 'w', 'MarkerSize', 5)
        PAx = polaraxes(Fig, 'Position', [0.75 (cnt-1)/nplt, 1/nplt 1/nplt]+[0.008 0 0 0]);
        polarplot(PAx, [inspect.angle.m(i), inspect.angle.m(i)], [0 1], '-', 'Color', standard_colors('red'), 'LineWidth', 3)
        PAx.RTick = [];
        PAx.ThetaTick = [];
        if cnt == nplt
            break
        end
    end
    Ax.XLim(2) = 100;
    exportgraphics(Fig, sprintf('./inspect/arophaseangle/autoreg_newfilt_fig%i.png', f), 'Resolution', 300)
    
end

%%

close all

session = {'placebo'};
aro_type = {'arousalemg'};
resolution = 20; % degrees
BinEdges = (-pi:resolution*pi/180:pi) - resolution*pi/360;

polar_resolution = 20; % degrees
PolBinEdges = (-pi:polar_resolution*pi/180:pi) - polar_resolution*pi/360;

for delay_s = [1 10]

    idx_sel_aw =  T.is_awakening & ismember(T.aro_type, aro_type) & ismember(T.ses, session);
    idx_sel_cs = ~T.is_awakening & ismember(T.aro_type, aro_type) & ismember(T.ses, session);

    Fig = figure();
    Fig.Color = 'w';
    Fig.Units = 'centimeters';
    Fig.Position(3:4) = [9, 7];
    Fig.Name = sprintf('-%i s', delay_s);
    % Create axes object
    Ax = axes('NextPlot', 'add');
    % Extract the phase angle of arousals leading to continued
    % sleep and awakenings and adjust phase angle for delay
    % - Continued sleep
    AData_cs = T.(sprintf('phase_d%i', delay_s))(idx_sel_cs);
    AData_cs = AData_cs + (delay_s/50)*2*pi; % Adjust for delay i.e., ISF phase angle + delay radians
    AData_cs(AData_cs > max(BinEdges)) = AData_cs(AData_cs > max(BinEdges))-2*pi;
    Pr_cs = histcounts(AData_cs, 'Normalization', 'probability', 'BinEdges', BinEdges);
    % - Awakenings
    AData_aw = T.(sprintf('phase_d%i', delay_s))(idx_sel_aw);
    AData_aw = AData_aw + (delay_s/50)*2*pi; % Adjust for delay i.e., ISF phase angle + delay radians
    AData_aw(AData_aw > max(BinEdges)) = AData_aw(AData_aw > max(BinEdges))-2*pi;
    Pr_aw = histcounts(AData_aw, 'Normalization', 'probability', 'BinEdges', BinEdges);
    % Plot ISF phase angle
    plot(Ax, -3*pi:pi/72:3*pi, 0.6*max([Pr_cs, Pr_aw])+cos((-3*pi:pi/72:3*pi))*0.4*max([Pr_cs, Pr_aw]),  '-', 'Color', [0.85 0.86 0.88], 'LineWidth', 5)
    % Plot bar graphs and set their face-color
    XData = BinEdges(1:end-1)+mean(diff(BinEdges)./2);
    XData = [XData-2*pi, XData, XData+2*pi]; %#ok<AGROW>
    YData = [Pr_cs', Pr_aw'];
    YData = [YData; YData; YData]; %#ok<AGROW>
    bh = bar(XData, YData, 'BarWidth', 1, 'LineStyle', 'none');
    bh(1).FaceColor = standard_colors('blue');
    bh(2).FaceColor = standard_colors('red');
    % Set axis properties
    Ax.Box = 'on';
    Ax.FontSize = 8;
    Ax.Layer = 'top';
    Ax.XLim = [0, 2*pi] + [-0.1*pi 0.1*pi];
    Ax.YLim = [0 2*max([Pr_cs, Pr_aw])];
    Ax.YTick = 0:0.02:Ax.YLim(2);
    Ax.XTick = (-3*pi:0.5*pi:3*pi);
    Ax.Color = [0.95 0.96 0.98];
    Ax.XGrid = 'on';
    Ax.TickLength = [0 0];
    Ax.XTickLabel = {'trough', 'asc', 'peak', 'desc', 'trough', 'asc', 'peak', 'desc', 'trough', 'asc', 'peak', 'desc', 'trough'};
    Ax.XLabel.String = 'Sigma ISF Phase Angle';
    Ax.XLabel.FontSize = 8;
    Ax.XLabel.FontWeight = 'bold';
    Ax.YLabel.String = 'Pr(arousal)';
    Ax.YLabel.FontSize = 8;
    Ax.YLabel.FontWeight = 'bold';
    % Insert legend
    legend(bh, {sprintf('continued sleep (n=%i)', sum(idx_sel_cs)), sprintf('awakening (n=%i)', sum(idx_sel_aw))}, 'Box', 'off');
    % Plot polaraxes insert
    Ax = polaraxes(Fig, 'NextPlot', 'add');
    clear ph;
    ph(1) = polarhistogram(AData_cs, 'Normalization', 'probability', 'BinEdges', PolBinEdges, 'LineStyle', 'none');
    ph(2) = polarhistogram(AData_aw, 'Normalization', 'probability', 'BinEdges', PolBinEdges, 'LineStyle', 'none'); %#ok<NASGU> 
    % Polar axis properties
    Ax.RTick = [];
    Ax.ThetaZeroLocation = 'top';
    Ax.RTickLabel = {};
    Ax.FontSize = 8;
    Ax.ThetaTick = [0 90 180 270];
    Ax.ThetaTickLabel = {'peak', 'desc', 'trough', 'asc'};
    Ax.Position = [0.2 0.6 0.25 0.25];
end

%%


fld = 'filt';
session = {'placebo'};
aro_type = {'arousalemg'};

idx_sel_aw =  T.is_awakening & ismember(T.aro_type, aro_type) & ismember(T.ses, session);
idx_sel_cs = ~T.is_awakening & ismember(T.aro_type, aro_type) & ismember(T.ses, session);


Fig = figure();
Ax = axes(Fig);
XData = -100:01:15;
idx_t = inspect.(fld).x{1} >= XData(1) & inspect.(fld).x{1} < XData(end);
XData = inspect.(fld).x{1}(idx_t);

YData_cs = cellfun(@(x, y) mean(y(:, x >= XData(1)-0.01 & x <= XData(end)+0.01), 1), inspect.(fld).x(idx_sel_cs), inspect.(fld).y(idx_sel_cs), 'UniformOutput', false);
YData_cs = cat(1, YData_cs{:})';
% EData_cs = tinv(0.975, sum(idx_sel_cs)-1).*std(YData_cs)./sqrt(sum(idx_sel_cs));
% YData_cs = mean(YData_cs, 1);

YData_aw = cellfun(@(x, y) mean(y(:, x >= XData(1)-0.01 & x <= XData(end)+0.01), 1), inspect.(fld).x(idx_sel_aw), inspect.(fld).y(idx_sel_aw), 'UniformOutput', false);
YData_aw = cat(1, YData_aw{:})';
% EData_aw = tinv(0.975, sum(idx_sel_aw)-1).*std(YData_aw)./sqrt(sum(idx_sel_aw));
% YData_aw = mean(YData_aw, 1);

hold on
linepatch(Ax, XData, YData_cs, 'EdgeAlpha', 0.1, 'EdgeColor', standard_colors('blue'));
linepatch(Ax, XData, YData_aw, 'EdgeAlpha', 0.1, 'EdgeColor', standard_colors('red'));
plot(Ax, XData, mean(YData_cs, 2), 'Color', standard_colors('blue'));
plot(Ax, XData, mean(YData_aw, 2), 'Color', standard_colors('red'));
plot(Ax, [0 0], [-8 8], '-w')
% patch(Ax, 'XData', [XData, fliplr(XData)], 'YData', [(YData_cs + EData_cs), fliplr(YData_cs - EData_cs)], ...
%     'LineStyle', 'none', ...
%     'FaceAlpha', 0.15, ...
%     'FaceColor', standard_colors('blue'))
% patch(Ax, 'XData', [XData, fliplr(XData)], 'YData', [(YData_aw + EData_aw), fliplr(YData_aw - EData_aw)], ...
%     'LineStyle', 'none', ...
%     'FaceAlpha', 0.15, ...
%     'FaceColor', standard_colors('red'))
% 
% plot(Ax, XData, YData_aw, 'Color', standard_colors('red'));

%% NON-AVERAGED ANGLES

%close all

this_delay = 0;
session = {'placebo'};
aro_type = {'arousalemg', 'arousal'};
resolution = 20; % degrees
BinEdges = (-pi:resolution*pi/180:pi) - resolution*pi/360;

polar_resolution = 20; % degrees
PolBinEdges = (-pi:polar_resolution*pi/180:pi) - polar_resolution*pi/360;
idx_sel_aw =  T.is_awakening & ismember(T.aro_type, aro_type) & ismember(T.ses, session);
idx_sel_cs = ~T.is_awakening & ismember(T.aro_type, aro_type) & ismember(T.ses, session);

Fig = figure();
Fig.Color = 'w';
Fig.Units = 'centimeters';
Fig.Position(3:4) = [9, 7];

% Create axes object
Ax = axes('NextPlot', 'add');
% Extract the phase angle of arousals leading to continued
% sleep and awakenings and adjust phase angle for delay
% - Continued sleep
AData_cs = T.phase_cell(idx_sel_cs);
AData_cs = [AData_cs{:}];
AData_cs = AData_cs + (this_delay/50)*2*pi; % Adjust for delay i.e., ISF phase angle + delay radians
AData_cs(AData_cs > max(BinEdges)) = AData_cs(AData_cs > max(BinEdges))-2*pi;
Pr_cs = histcounts(AData_cs, 'Normalization', 'probability', 'BinEdges', BinEdges);
% - Awakenings
AData_aw = T.phase_cell(idx_sel_aw);
AData_aw = [AData_aw{:}];
AData_aw = AData_aw + (this_delay/50)*2*pi; % Adjust for delay i.e., ISF phase angle + delay radians
AData_aw(AData_aw > max(BinEdges)) = AData_aw(AData_aw > max(BinEdges))-2*pi;
Pr_aw = histcounts(AData_aw, 'Normalization', 'probability', 'BinEdges', BinEdges);
% Plot ISF phase angle
plot(Ax, -3*pi:pi/72:3*pi, 0.6*max([Pr_cs, Pr_aw])+cos((-3*pi:pi/72:3*pi))*0.4*max([Pr_cs, Pr_aw]),  '-', 'Color', [0.85 0.86 0.88], 'LineWidth', 5)
% Plot bar graphs and set their face-color
XData = BinEdges(1:end-1)+mean(diff(BinEdges)./2);
XData = [XData-2*pi, XData, XData+2*pi]; 
YData = [Pr_cs', Pr_aw'];
YData = [YData; YData; YData];
bh = bar(XData, YData, 'BarWidth', 1, 'LineStyle', 'none');
bh(1).FaceColor = standard_colors('blue');
bh(2).FaceColor = standard_colors('red');
% Set axis properties
Ax.Box = 'on';
Ax.FontSize = 8;
Ax.Layer = 'top';
Ax.XLim = [0, 2*pi] + [-0.08*pi 0.08*pi];
Ax.YLim = [0 2*max([Pr_cs, Pr_aw])];
Ax.YTick = 0:0.02:Ax.YLim(2);
Ax.XTick = (-3*pi:0.5*pi:3*pi);
Ax.Color = [0.95 0.96 0.98];
Ax.XGrid = 'on';
Ax.TickLength = [0 0];
Ax.XTickLabel = {'trough', 'asc', 'peak', 'desc', 'trough', 'asc', 'peak', 'desc', 'trough', 'asc', 'peak', 'desc', 'trough'};
Ax.XLabel.String = 'Sigma ISF Phase Angle';
Ax.XLabel.FontSize = 8;
Ax.XLabel.FontWeight = 'bold';
Ax.YLabel.String = 'Pr(arousal)';
Ax.YLabel.FontSize = 8;
Ax.YLabel.FontWeight = 'bold';
% Insert legend
legend(bh, {sprintf('continued sleep (n=%i)', sum(idx_sel_cs)), sprintf('awakening (n=%i)', sum(idx_sel_aw))}, 'Box', 'off');
% Plot polaraxes insert
Ax = polaraxes(Fig, 'NextPlot', 'add');
clear ph;
ph(1) = polarhistogram(AData_cs, 'Normalization', 'probability', 'BinEdges', PolBinEdges, 'LineStyle', 'none');
ph(2) = polarhistogram(AData_aw, 'Normalization', 'probability', 'BinEdges', PolBinEdges, 'LineStyle', 'none');
% Polar axis properties
Ax.RTick = [];
Ax.ThetaZeroLocation = 'top';
Ax.RTickLabel = {};
Ax.FontSize = 8;
Ax.ThetaTick = [0 90 180 270];
Ax.ThetaTickLabel = {'peak', 'desc', 'trough', 'asc'};
Ax.Position = [0.2 0.6 0.25 0.25];

fname = struct();
fname.fcut = strrep(sprintf('%.4f%.4f',filtcfg.cutoff), '0.', 'x');
fname.tbw = sprintf('1x%i', 1/filtcfg.transbw);
fname.rdev = strrep(sprintf('%.4f',filtcfg.rippledev), '0.', 'x');
exportgraphics(Fig, sprintf('./inspect/isfphase-nonav_ftype-%s_forder-%i_fcut-%s_tbw-%s_rdev-%s_atype-%s.png', filtcfg.wintype, filtcfg.order, fname.fcut, fname.tbw, fname.rdev, strrep(append_type, '_', '')), 'Resolution', 300)
