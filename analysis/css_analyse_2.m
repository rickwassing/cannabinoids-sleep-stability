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


%% Define the set of parietal channels to focus on
chanlocs = template_to_chanlocs(which('GSN-HydroCel-257.sfp'));
[~, system] = ishdeeg({chanlocs.labels});
incl = hdeeg_scalpchannels(system);
chanlocs = chanlocs(ismember({chanlocs.labels}, incl));
chanlocs = channel_clusters(chanlocs, 'mff');

YData = arrayfun(@(s) s.features(2).data, ISF, 'UniformOutput', false);
YData = cat(2, YData{:});
YData = mean(YData, 2);

idx_chan = [77, 78, 79, 85, 86, 87, 92, 93, 94, 99, 104, 110, 111, 112, 113, 120, 121, 122];
chanlabels = {chanlocs(idx_chan).labels};

Fig = figure();
Fig.Units = 'centimeters';
Fig.Position(3:4) = [8 8];
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Plot as topoplot
topoplot(YData, chanlocs, ...
    'colormap', batlow, ...
    'hlinewidth', 1, ...
    'maplimits', [0, 2], ...
    'numcontour', 0, ...
    'conv', 'on', ...
    'whitebk', 'on', ...
    'emarker2', {idx_chan, '.', 'w', 15, 1});

exportgraphics(Fig, 'figures/supp_selparietalchans.png', 'Resolution', 600)

%%
Files = dir('derivatives/EEG-segmented/sub-*/ses-*/sub-*_desc-sigmanrembout*.set');

clear SIGMA
for i = 1:length(Files)
    SIGMA(i) = LoadDataset(fullfile(Files(i).folder, Files(i).name), 'all');
    SIGMA(i) = pop_select(SIGMA(i), 'channel', idx_chan);
end

%% Emperical values for the filter cutoffs
Files = dir('derivatives/EEG-output-fstlvl/sub-*/ses-*/sub-*_desc-a1cabssigma_fstlvl.mat');
fcutoff = [];
for i = 1:length(Files)
    tmp = LoadDataset(fullfile(Files(i).folder, Files(i).name), 'matrix');
    bw = [...
        tmp.features(1).data(idx_chan) - tmp.features(3).data(idx_chan)/2, ...
        tmp.features(1).data(idx_chan) + tmp.features(3).data(idx_chan)/2];
    fcutoff = [fcutoff; mean(bw)]; %#ok<AGROW>
end
fcutoff = mean(fcutoff);
%%
H = struct();
% For each bout, filter the full bout, and filter the bout minus the filter order

ANGA = [];
ANGB = [];
DANG = [];
DAMP = [];

for i = 1:length(SIGMA)

    bouts = find(strcmpi({SIGMA(i).event.type}, 'boundary'));
    bouts = ceil([1, [SIGMA(i).event(bouts).latency], SIGMA(i).pnts]);

    forder = pop_firwsord('hamming', SIGMA(i).srate, fcutoff(1));
    if mod(forder, 2) == 1
        forder = forder+1;
    end
    
    % Manually overwrite filter order
    forder = 1000;

    for b = 1:length(bouts)-1
        if bouts(b+1)-bouts(b) < 6*forder
            continue
        end

        tmp_a = pop_select(SIGMA(i), 'point', [bouts(b), bouts(b+1)]);
        tmp_a.times = tmp_a.times./1000;
        tmp_a.data = detrend(tmp_a.data', 0, 'omitnan')';
        tmp_a = pop_firws(tmp_a, ...
            'fcutoff', fcutoff, ...
            'ftype', 'bandpass', ...
            'wtype', 'hamming', ...
            'forder', forder, ...
            'minphase', 0);

        tmp_b = pop_select(SIGMA(i), 'point', [bouts(b), bouts(b+1)]+[0 -2*forder]);
        tmp_b.times = tmp_b.times./1000;
        tmp_b.data = detrend(tmp_b.data', 0, 'omitnan')';
        tmp_b.append = flip(tmp_b.data(:, end-2*forder+1:end), 2);
        if true
            tmp_b.data = [tmp_b.data, tmp_b.append];
            tmp_b.pnts = size(tmp_b.data, 2);
            tmp_b.xmin = 0;
            tmp_b.xmax = tmp_b.pnts/tmp_b.srate;
            tmp_b.times = linspace(tmp_b.xmin, tmp_b.xmax, tmp_b.pnts);
        end
        tmp_b = pop_firws(tmp_b, ...
            'fcutoff', fcutoff, ...
            'ftype', 'bandpass', ...
            'wtype', 'hamming', ...
            'forder', forder, ...
            'minphase', 0);
        if true
            tmp_b.data = tmp_b.data(:, 1:end-2*forder); % cut away the appended data
        end

        H.a.sig = mean(tmp_a.data, 1)';
        H.a.sig = detrend(H.a.sig, 0);
        H.a.x = hilbert(H.a.sig);
        H.a.sig = H.a.sig(1:end-2*forder);
        H.a.ang = angle(H.a.x(1:end-2*forder));
        H.a.amp = abs(H.a.x(1:end-2*forder));

        H.b.sig = mean(tmp_b.data, 1)';
        H.b.sig = detrend(H.b.sig, 0);
        H.b.x = hilbert(H.b.sig);
        H.b.ang = angle(H.b.x);
        H.b.amp = abs(H.b.x);

        dang = asrow(circ_dist(H.a.ang, H.b.ang));
        damp = asrow(abs(H.a.amp - H.b.amp));

        % store the first 'forder' samples, and the last
        ANGA = [ANGA; asrow(H.a.sig(end-4*forder+1:end))];
        ANGB = [ANGB; asrow(H.b.sig(end-4*forder+1:end))];
        DANG = [DANG; dang(end-4*forder+1:end)];
        DAMP = [DAMP; damp(end-4*forder+1:end)];

    end
end

%%
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
Ax(ai) = axes();
XData = 1:2;
YData = [4000, 0, 0, 0; 0, 0, 4000, 2000];
b = barh(XData, YData, 'stacked', 'BarWidth', 0.5);
b(1).FaceColor = standard_colors('red');
b(2).FaceColor = standard_colors('red').^0.2;
b(3).FaceColor = standard_colors('blue');
b(4).FaceColor = standard_colors('blue').^0.2;

% Axis props
Ax(ai).Box = 'off';
Ax(ai).OuterPosition = [0/12 1-4/12 12/12 4/12] + [0 0 0.055 0];
Ax(ai).YTickLabel = {'test signal', 'control signal'};
Ax(ai).TickLength = [0, 0];
Ax(ai).XTick = [0, 4000, 6000];
Ax(ai).XTickLabel = {'-400', '0', '+200'};
Ax(ai).XLabel.String = 'time (s)';
Ax(ai).XAxisLocation = 'top';

% -------------------------------------------------------------------------
% Plot the control and test signals
ai = ai+1;
Ax(ai) = axes('NextPlot', 'add', 'Layer', 'top');

patch(Ax(ai), 'XData', [3500 4000 4000 3500], 'YData', [-1 -1 1 1], 'LineStyle', 'none', 'FaceColor', [0.92 0.93 0.95])
plot(1:size(ANGA,2), mean(ANGA, 1), '-', 'LineWidth', 1, 'Color', standard_colors('blue'))
plot(1:size(ANGA,2), mean(ANGB, 1), '-', 'LineWidth', 1, 'Color', standard_colors('red'))

% Axis props
Ax(ai).Box = 'on';
Ax(ai).OuterPosition = [0/12 1-7/12 8/12 3/12]+[0.0425 0 0.045 0];
Ax(ai).TickLength = [0, 0];
Ax(ai).XTick = [0, 3500, 4000];
Ax(ai).XTickLabel = {'-400', '-50', '0'};
Ax(ai).YLim = [-0.2 0.2];
Ax(ai).YTick = [];

% -------------------------------------------------------------------------
% Plot a Zoom-in of the control and test signals
ai = ai+1;
Ax(ai) = axes('NextPlot', 'add', 'Color', [0.92 0.93 0.95]);

plot(mean(ANGA(:, end-500:end), 1), '-', 'LineWidth', 1, 'Color', standard_colors('blue'))
plot(mean(ANGB(:, end-500:end), 1), '-', 'LineWidth', 1, 'Color', standard_colors('red'))

% Axis props
Ax(ai).Box = 'on';
Ax(ai).OuterPosition = [8/12 1-7/12 4/12 3/12]+[0.025 0 -0.04 0];
Ax(ai).TickLength = [0, 0];
Ax(ai).XLim = [0, 500];
Ax(ai).XTick = 0:50:500;
Ax(ai).XTickLabel = {'-50', '', '-40', '', '-30', '', '-20', '', '-10', '', '0'};
Ax(ai).XTickLabelRotation = 0;
Ax(ai).XGrid = 'on';
Ax(ai).YLim = Ax(ai-1).YLim;
Ax(ai).YTick = [];

% -------------------------------------------------------------------------
% Plot the difference in angle
ai = ai+1;
Ax(ai) = axes('NextPlot', 'add', 'Color', [0.92 0.93 0.95]);

YData = DANG(:, end-500:end)';
XData = repmat([0:500, nan], 1, size(YData, 1));
YData = [YData; nan(1, size(YData, 2))];
YData = YData(:);

patch('XData', XData, 'YData', YData, 'EdgeColor', 'k', 'EdgeAlpha', 0.08);

% Axis props
Ax(ai).Box = 'on';
Ax(ai).OuterPosition = [8/12 1-12/12 4/12 5/12]+[0.025 0 -0.04 0];
Ax(ai).TickLength = [0, 0];
Ax(ai).XLim = [0, 500];
Ax(ai).XTick = [0, 500];
Ax(ai).XTickLabel = {'-50', '0'};
Ax(ai).YTick = [-pi, 0, pi];
Ax(ai).YTickLabel = {'-\pi', '0', '\pi'};
Ax(ai).YLim = [-3.6, 3.6];


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
    ai = ai + 1;
    Ax(ai) = polaraxes('NextPlot', 'add', 'Color', [0.92 0.93 0.95]);
    Ax(ai).OuterPosition = [(i)*1.175/12 1-12/12 1.2/12 5/12] + [0.01 0 0 0];
    Ax(ai).ThetaTick = [0, 90, 180, 270];
    Ax(ai).ThetaTickLabel = {};
    Ax(ai).RTick = [];
    Ax(ai).Title.String = sprintf('-%i to -%i s', delay(i)/10, delay(i+1)/10);
    Ax(ai).Title.FontSize = 8;
    Ax(ai).Title.FontWeight = 'normal';
    Theta = DANG(:, (end-delay(i)+1):(end-delay(i+1)));
    Theta = Theta(:);
    h = polarhistogram(Ax(ai), Theta, 36, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', [0.5 0.5 0.5]);
    polarplot(Ax(ai), [circ_mean(Theta), circ_mean(Theta)], [0 max(h.Values)], '-k', 'LineWidth', 1.5)
    Ax(ai).RLim = [0, max(h.Values)];
end

exportgraphics(Fig, 'figures/supp_2.png', 'Resolution', 300)
disp('Done saving')
%%



ai = ai+1;
Ax(ai) = axes();
Ax = axes('NextPlot', 'add');

YData = [DANG'; nan(1, size(DANG, 1))];
YData = YData(:);
XData = repmat([1:4000, nan], 1, size(DANG,1));
patch('XData', XData, 'YData', YData, 'EdgeColor', 'k', 'EdgeAlpha', 0.08);
plot(circ_mean(DANG, [], 1), '-k', 'LineWidth', 1.5)
plot(circ_mean(DANG, [], 1)+circ_std(DANG, [], [], 1), ':k', 'LineWidth', 1.5)
plot(circ_mean(DANG, [], 1)-circ_std(DANG, [], [], 1), ':k', 'LineWidth', 1.5)
plot([0, 4*forder], [-pi/8, -pi/8], '-r')
plot([0, 4*forder], [pi/8, pi/8], '-r')
Ax.XTick = 1:50:4*forder;
Ax.XTickLabel = (Ax.XTick-4*forder-1)./10;
Ax.YTick = -pi:0.25*pi:pi;


% polarhistogram(circ_mean(DANG(:, end-100:end-50),[], 2), 36)
