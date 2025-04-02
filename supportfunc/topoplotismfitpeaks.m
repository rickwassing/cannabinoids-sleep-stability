addpath('/Users/rickwassing/Local/eeglab/latest'); eeglab; close all;
addpath('/Users/rickwassing/Local/fieldtrip/latest'); ft_defaults;
addpath(genpath('/Users/rickwassing/Local/EEG_Processor/develop'));

%%

CMap = load('colormap_roma.mat');
CMap = CMap.roma;

%%

Files = dir('derivatives/EEG-preproc/sub-*/ses-*/sub-*_desc-ismfit_struct.mat');

%%

fname = 'derivatives/EEG-preproc/sub-r0021jb/ses-1/sub-r0021jb_ses-1_task-psg_run-1_desc-nrem_channels.tsv';
specchans = readSidecarTSV(fname, 'channels');
specchans(~strcmpi(specchans.type, 'EEG'), :) = [];
chanlocs = readlocs('GSN-HydroCel-257.sfp');
chanlocs(~ismember({chanlocs.labels}, specchans.name)) = [];

%%
FIT = struct();
for i = 1:length(Files)
    tmp = load(fullfile(Files(i).folder, Files(i).name));
    FIT(i).f = tmp.FIT;
end

%%

for i = 1:length(Files)
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % YData
    fld = 'mu';
    YData = [FIT(i).f.(fld)];
    YData(end) = []; % Remove average
    switch fld
        case 'peak'
            if any(YData > 4)
                MaxYData = min(YData(YData > 4));
                YData(YData > 4) = 4;
            else
                MaxYData = 4;
            end
            contourvals = (YData > 1) + (YData > 2) + (YData > 3);
            CLim = [0, 4];
            Ticks = sort([0, 1, 2, 3, 4]);
            TickLabels = [0, 1, 2, 3, MaxYData];
            Label = 'peak amplitude (norm.)';
            outfilepath = 'group-level/topopeaks';
        case 'mu'
            if any(YData > 0.04)
                MaxYData = min(YData(YData > 0.04));
                YData(YData > 0.04) = 0.04;
            else
                MaxYData = 0.04;
            end
            contourvals = (YData > 0.01) + (YData > 0.02) + (YData > 0.03);
            CLim = [0, 0.04];
            Ticks = sort([0, 0.01, 0.02, 0.03, 0.04]);
            TickLabels = [0, 0.01, 0.02, 0.03, MaxYData];
            Label = 'peak frequency (Hz)';
            outfilepath = 'group-level/topofreqs';
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Create figure
    Fig = figure('Color', 'w', 'Position', [560, 640, 290 200]);
    Ax = axes();
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Draw topoplot

    topoplot(YData, chanlocs, ...
        'headrad', 0.575, ...
        'whitebk', 'on', ...
        'conv', 'on', ...
        'numcontour', length(unique(contourvals))-1, ...
        'contourvals', contourvals);
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Reset figure and axis properties
    Fig.Color = 'w';
    Fig.Colormap = CMap;
    Ax.CLim = CLim;
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Colobar
    CBar = colorbar();
    CBar.Ticks = Ticks;
    CBar.TickLabels = TickLabels;
    CBar.Label.String = Label;
    CBar.FontSize = 8;
    CBar.Label.FontSize = 10;
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Export graphics
    outfilename = strrep(Files(i).name, 'desc-ismfit_struct.mat', 'ismpeak.png');
    if exist(outfilepath, 'dir') == 0
        mkdir(outfilepath);
    end
    exportgraphics(Fig, fullfile(outfilepath, outfilename));

end


%%

fld = 'peak';
YData = arrayfun(@(f) [f.f.(fld)], FIT, 'UniformOutput', false);
YData = cat(1, YData{:});
YData = mean(YData);

switch fld
    case 'peak'
        if any(YData > 4)
            MaxYData = min(YData(YData > 4));
            YData(YData > 4) = 4;
        else
            MaxYData = 4;
        end
        contourvals = (YData > 1) + (YData > 2) + (YData > 3);
        CLim = [0.5, 3.5];
        Ticks = sort([0, 1, 2, 3, 4]);
        TickLabels = [0, 1, 2, 3, MaxYData];
        Label = 'peak amplitude (norm.)';
        outfilepath = 'group-level/topopeaks';
    case 'mu'
        if any(YData > 0.04)
            MaxYData = min(YData(YData > 0.04));
            YData(YData > 0.04) = 0.04;
        else
            MaxYData = 0.04;
        end
        contourvals = (YData > 0.01) + (YData > 0.02) + (YData > 0.03);
        CLim = [0.01, 0.03];
        Ticks = sort([0, 0.01, 0.02, 0.03, 0.04]);
        TickLabels = [0, 0.01, 0.02, 0.03, MaxYData];
        Label = 'peak frequency (Hz)';
        outfilepath = 'group-level/topofreqs';
end
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Create figure
Fig = figure('Color', 'w', 'Position', [560, 640, 290 200]);
Ax = axes();
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Draw topoplot
topoplot(YData, chanlocs, ...
    'headrad', 0.575, ...
    'whitebk', 'on', ...
    'conv', 'on', ...
    'numcontour', length(unique(contourvals))-1, ...
    'contourvals', contourvals);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Reset figure and axis properties
Fig.Color = 'w';
Fig.Colormap = CMap;
Ax.CLim = CLim;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Colobar
CBar = colorbar();
CBar.Ticks = Ticks;
CBar.TickLabels = TickLabels;
CBar.Label.String = Label;
CBar.FontSize = 8;
CBar.Label.FontSize = 10;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Export graphics
outfilename = '01_average_peak.png';
if exist(outfilepath, 'dir') == 0
    mkdir(outfilepath);
end
exportgraphics(Fig, fullfile(outfilepath, outfilename));
%%

sigma_filepath = 'derivatives/EEG-preproc/sub-r0091/ses-1/sub-r0091_ses-1_task-psg_run-1_desc-sigma_eeg.set';
fit_filepath = 'derivatives/EEG-preproc/sub-r0091/ses-1/sub-r0091_ses-1_task-psg_run-1_desc-ismfit_struct.mat';
SIGMA = LoadDataset(sigma_filepath, 'all');
hr = SIGMA.data(end, :);
filtorder = pop_firwsord('hamming', SIGMA.srate, 0.01);
INFRA = pop_firws(SIGMA, ...
    'fcutoff', [0.005, 0.055], ...
    'ftype', 'bandpass', ...
    'wtype', 'hamming', ...
    'forder', filtorder, ...
    'minphase', 0);
idx_nan = sum(isnan(INFRA.data), 1) > 0;
INFRA.data(:, idx_nan) = 0;
FIT = load(fit_filepath);
FIT = FIT.FIT;
idx_chan = [FIT.peak] > 1;
idx_chan(end) = []; % remove the last element (the average)
infra_mu = mean(INFRA.data(idx_chan, :), 1);
tmp = hilbert(infra_mu);
amp = abs(tmp)';
phs = angle(tmp)';
infra_mu(idx_nan) = NaN;
amp(idx_nan) = NaN;
phs(idx_nan) = NaN;

%%

Fig = figure('Color', 'w')';
Ax = axes(...
    'Box', 'on', ...
    'NextPlot', 'add', ...
    'Color', [233, 236, 239]./255, ...
    'XColor', [0.5, 0.5, 0.5], ...
    'YColor', [0.5, 0.5, 0.5], ...
    'XGrid', 'on', ...
    'YGrid', 'off', ...
    'YMinorGrid', 'off', ...
    'GridAlpha', 0.25, ...
    'LineWidth', 0.5, ...
    'TickLength', [0, 0], ...
    'FontSize', 10, ...
    'Layer', 'top', ...
    'Units', 'normalized');

hold on
plot(SIGMA.times, nanzscore(mean(SIGMA.data(idx_chan, :), 1, 'omitnan')))
plot(SIGMA.times, nanzscore(infra_mu)+2, '-', 'LineWidth', 1.5)

for i = 1:length(SIGMA.event)
    if ~strcmpi(SIGMA.event(i).type, 'arousalemg')
        continue
    end
    XData = [SIGMA.event(i).latency, SIGMA.event(i).latency+SIGMA.event(i).duration] ./ SIGMA.srate;
    YData = [7, 7];
    plot(XData, YData, '-', 'Color', CMap(end, :), 'LineWidth', 4)
end
