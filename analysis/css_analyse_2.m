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
        H.a.ang = angle(H.a.x(1:end-2*forder));
        H.a.amp = abs(H.a.x(1:end-2*forder));

        H.b.sig = mean(tmp_b.data, 1)';
        H.b.sig = detrend(H.b.sig, 0);
        H.b.x = hilbert(H.b.sig);
        H.b.ang = angle(H.b.x);
        H.b.amp = abs(H.b.x);

        dang = asrow(mod(H.a.ang - H.b.ang + pi, 2*pi) - pi);
        damp = asrow(abs(H.a.amp - H.b.amp));

        % store the first 'forder' samples, and the last
        DANG = [DANG; dang(end-4*forder+1:end)];
        DAMP = [DAMP; damp(end-4*forder+1:end)];

    end
end

%%
check = double(DANG(:, end));

%%
close all
figure

YData = [DANG'; nan(1, size(DANG, 1))];
YData = YData(:);
XData = repmat([1:4000, nan], 1, size(DANG,1));
Ax = axes('NextPlot', 'add');
patch('XData', XData, 'YData', YData, 'EdgeColor', 'k', 'EdgeAlpha', 0.08);
plot(mean(DANG), '-k', 'LineWidth', 1.5)
plot(prctile(DANG, 95), ':k', 'LineWidth', 1.5)
plot(prctile(DANG, 5), ':k', 'LineWidth', 1.5)
plot([0, 4*forder], [-pi/8, -pi/8], '-r')
plot([0, 4*forder], [pi/8, pi/8], '-r')
Ax.XTick = 1:50:4*forder;
Ax.XTickLabel = (Ax.XTick-4*forder-1)./10;
Ax.YTick = -pi:0.25*pi:pi;

figure;

polarhistogram(circ_mean(DANG(:, end-100:end-50),[], 2), 36)
