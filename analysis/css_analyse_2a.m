% function css_analyse_2a
clc
clear
close all
% -------------------------------------------------------------------------
% Determine the impact of the filter edge artefact on the phase estimate of
% the ISF
% -------------------------------------------------------------------------
% Define the set of parietal channels to focus on
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
% Load the SIGMA power timeseries
Files = dir('derivatives/EEG-segmented/sub-*/ses-*/sub-*_desc-sigmanrembout*.set');
clear SIGMA
for i = 1:length(Files)
    SIGMA(i) = LoadDataset(fullfile(Files(i).folder, Files(i).name), 'all'); %#ok<SAGROW> 
    SIGMA(i) = pop_select(SIGMA(i), 'channel', idx_chan); %#ok<SAGROW> % Select parietal channels
end
% -------------------------------------------------------------------------
% Get emperical values for the filter cutoffs
filtcfg.cutoff = [0.0138 0.0259];
filtcfg.wintype = 'kaiser';
filtcfg.transbw = 1/50;
filtcfg.rippledev = 0.05;
filtcfg.warg = 2;
filtcfg.adj = [0 0];%[0.005 -0.001];
filtcfg.order = pop_firwsord(filtcfg.wintype, SIGMA(1).srate, filtcfg.transbw, filtcfg.rippledev); % transition bandwidth of 0.01 Hz
filtcfg.order
% -------------------------------------------------------------------------
% Create two copies of the sigma-timeseries. One which is 6 times the
% filter length, and another that is 4 times the filter length (cropped
% from the end). Then filter these timeseries, crop the long timeseries to
% the same size as the shorter one, and apply Hilbert transform to get the
% instantaneous phase estimate. The difference in the phase estimate is the
% impact of the edge artifact of the filter.
% -------------------------------------------------------------------------
% Init
H = struct();
ANGA = [];
ANGB = [];
SIGA = [];
SIGB = [];
DANG = [];
DAMP = [];
% -------------------------------------------------------------------------
% Loop over files and extract bouts
cbouts = struct();
cbouts.total = 0;
cbouts.tooshort = 0;
cbouts.analysed = 0;
trem = now();
for i = 1:length(SIGMA)
    % ---------------------------------------------------------------------
    % Get indices of the selected NREM bouts
    bouts = find(strcmpi({SIGMA(i).event.type}, 'boundary'));
    bouts = ceil([1, [SIGMA(i).event(bouts).latency], SIGMA(i).pnts]);
    % ---------------------------------------------------------------------
    % Loop over bouts...
    for b = 1:length(bouts)-1
        % ... continue to the next bout if is too short
        cbouts.total = cbouts.total+1;
        if bouts(b+1)-bouts(b) < 6*filtcfg.order
            cbouts.tooshort = cbouts.tooshort+1;
            continue
        end
        cbouts.analysed = cbouts.analysed+1;
        % -----------------------------------------------------------------
        % Create control signal
        sig_control = pop_select(SIGMA(i), 'point', [bouts(b), bouts(b+1)]);
        sig_control.times = sig_control.times./1000;
        sig_control.data = detrend(sig_control.data', 0, 'omitnan')'; % demean
        % -----------------------------------------------------------------
        % Create test signal
        sig_test = pop_select(SIGMA(i), 'point', [bouts(b), bouts(b+1)]+[0 -ford_mult*filtcfg.order]);
        sig_test.times = sig_test.times./1000;
        sig_test.data = detrend(sig_test.data', 0, 'omitnan')'; % demean
        % -----------------------------------------------------------------
        % Get prepend and append signals
        sig_test = signalappend(sig_test, append_type, ford_mult, filtcfg);
        sig_control = signalappend(sig_control, append_type, ford_mult, filtcfg);
        % -----------------------------------------------------------------
        % Prepend and append data to cover for filter edge artefact
        sig_test.data = [sig_test.prepend, sig_test.data, sig_test.append];
        sig_test.pnts = size(sig_test.data, 2);
        sig_test.xmin = 0;
        sig_test.xmax = sig_test.pnts/sig_test.srate;
        sig_test.times = linspace(sig_test.xmin, sig_test.xmax, sig_test.pnts);
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        sig_control.data = [sig_control.prepend, sig_control.data, sig_control.append];
        sig_control.pnts = size(sig_control.data, 2);
        sig_control.xmin = 0;
        sig_control.xmax = sig_control.pnts/sig_control.srate;
        sig_control.times = linspace(sig_control.xmin, sig_control.xmax, sig_control.pnts);
        % -----------------------------------------------------------------
        % Filter the data
        sig_control = pop_firws(sig_control, ...
            'fcutoff', filtcfg.cutoff+filtcfg.adj, ...
            'ftype', 'bandpass', ...
            'wtype', filtcfg.wintype, ...
            'warg', filtcfg.warg, ...
            'forder', filtcfg.order, ...
            'plotfresp', false, ...
            'minphase', 0);
        sig_test = pop_firws(sig_test, ...
            'fcutoff', filtcfg.cutoff+filtcfg.adj, ...
            'ftype', 'bandpass', ...
            'wtype', filtcfg.wintype, ...
            'warg', filtcfg.warg, ...
            'forder', filtcfg.order, ...
            'plotfresp', false, ...
            'minphase', 0);
        % -----------------------------------------------------------------
        % cut away the appended data
        sig_control.data = sig_control.data(:, length(sig_control.prepend)+1:end-length(sig_control.append));
        sig_control.pnts = size(sig_control.data, 2);
        sig_control.xmin = 0;
        sig_control.xmax = sig_control.pnts/sig_control.srate;
        sig_control.times = linspace(sig_control.xmin, sig_control.xmax, sig_control.pnts);
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        sig_test.data = sig_test.data(:, length(sig_test.prepend)+1:end-length(sig_test.append));
        sig_test.pnts = size(sig_test.data, 2);
        sig_test.xmin = 0;
        sig_test.xmax = sig_test.pnts/sig_test.srate;
        sig_test.times = linspace(sig_test.xmin, sig_test.xmax, sig_test.pnts);
        % -----------------------------------------------------------------
        % Apply Hilbert to control signal
        H.control.sig = sig_control.data(:, 1:end-ford_mult*filtcfg.order)';
        H.control.sig = detrend(H.control.sig, 0);
        H.control.x = hilbert(H.control.sig);
        H.control.ang = angle(H.control.x);
        H.control.amp = abs(H.control.x);
        % Apply Hilbert to test signal
        H.test.sig = sig_test.data';
        H.test.sig = detrend(H.test.sig, 0);
        H.test.x = hilbert(H.test.sig);
        H.test.ang = angle(H.test.x);
        H.test.amp = abs(H.test.x);
        % Caclulate the difference in phase and amplitude
        dang = circ_dist(H.control.ang, H.test.ang);
        damp = abs(H.control.amp - H.test.amp);
        % Store estimated angle and difference in angle for this bout
        SIGA = [SIGA; H.control.sig(end-4*filtcfg.order+1:end, :)']; %#ok<AGROW>
        SIGB = [SIGB; H.test.sig(end-4*filtcfg.order+1:end, :)']; %#ok<AGROW>
        ANGA = [ANGA; H.control.ang(end-4*filtcfg.order+1:end, :)']; %#ok<AGROW>
        ANGB = [ANGB; H.test.ang(end-4*filtcfg.order+1:end, :)']; %#ok<AGROW>
        DANG = [DANG; dang(end-4*filtcfg.order+1:end, :)']; %#ok<AGROW>
        DAMP = [DAMP; damp(end-4*filtcfg.order+1:end, :)']; %#ok<AGROW>
    end
    trem = remainingTime(trem, length(SIGMA), true);
end
%% -------------------------------------------------------------------------
% Create supplementary figure
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
bh(2).FaceColor = standard_colors('red');
bh(3).FaceColor = standard_colors('red').^0.6;
bh(4).FaceColor = standard_colors('white');
bh(5).FaceColor = standard_colors('white');
bh(6).FaceColor = standard_colors('blue');
bh(7).FaceColor = standard_colors('blue').^0.6;
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
plot(1:size(SIGA, 2), SIGB(idx_sel, :).*4+5, '-', 'LineWidth', 1, 'Color', standard_colors('red'))
plot(1:size(SIGA, 2), SIGA(idx_sel, :).*4+5, '-', 'LineWidth', 1, 'Color', standard_colors('blue'))
plot(1:size(ANGB, 2), ANGB(idx_sel, :)-pi, '-', 'LineWidth', 1, 'Color', standard_colors('red'))
plot(1:size(ANGA, 2), ANGA(idx_sel, :)-pi, '-', 'LineWidth', 1, 'Color', standard_colors('blue'))
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
plot(1:size(SIGA, 2), SIGB(idx_sel, :).*4+5, '-', 'LineWidth', 1, 'Color', standard_colors('red'))
plot(1:size(SIGA, 2), SIGA(idx_sel, :).*4+5, '-', 'LineWidth', 1, 'Color', standard_colors('blue'))
plot(1:size(ANGB, 2), ANGB(idx_sel, :)-pi, '-', 'LineWidth', 1, 'Color', standard_colors('red'))
plot(1:size(ANGA, 2), ANGA(idx_sel, :)-pi, '-', 'LineWidth', 1, 'Color', standard_colors('blue'))
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
Ax(ai).OuterPosition = [0 1-12/12 1.2/12 5/12] + [0 -0.0125 0.025 0.025];
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
    h = polarhistogram(Ax(ai), Theta, 36, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', [0.5 0.5 0.5]);
    polarplot(Ax(ai), [circ_mean(Theta), circ_mean(Theta)], [0 max(h.Values)], '-k', 'LineWidth', 1.5)
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    Ax(ai).RLim = [0, max(h.Values)];
end
% -------------------------------------------------------------------------
% Save image
exportgraphics(Fig, sprintf('figures/supp_2_append-%s.png', append_type), 'Resolution', 300)
disp('Done saving')

% end