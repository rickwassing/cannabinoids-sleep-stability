cd('/Volumes/sleep/Sleep/3. ACTIVE STUDIES/CUPID/Arousal paper backup_CANSLEEP')
clc
PHEN = readtable('Scoring log and notes_CANSLEEP arousal.xlsx');

Files = dir('derivatives/EEG-preproc/sub-*/ses-*/sub*ismfit*');

Fig = figure('Color', 'w', 'Position', [560, 640, 290 200]);
Ax = axes('NextPlot', 'add', 'FontSize', 10, 'XTick', 0:0.01:0.1, 'YTick', 1:10, 'XGrid', 'on', 'TickLength', [0 0], 'Box', 'on', 'Color', [0.96 0.97 0.99]);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
H.data = patch(Ax, 'XData', NaN, 'YData', NaN, 'EdgeColor', 'k', 'EdgeAlpha', 0.1, 'FaceColor', 'none');
H.peak = plot(Ax, NaN, NaN, '.', 'Color', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 6);
H.avdata = plot(Ax, NaN, NaN, '-w', 'LineWidth', 1.5);

% H.gaus(1) = plot(Ax, NaN, NaN, ':r', 'LineWidth', 0.5);
% H.gaus(2) = plot(Ax, NaN, NaN, ':r', 'LineWidth', 0.5);
% H.gaus(3) = plot(Ax, NaN, NaN, ':r', 'LineWidth', 0.5);
% H.gaus(4) = plot(Ax, NaN, NaN, ':r', 'LineWidth', 0.5);
% H.thres = plot(Ax, NaN, NaN, '--k', 'LineWidth', 0.75);
% H.lolim = plot(Ax, NaN, NaN, '--r', 'LineWidth', 0.75);
% H.hilim = plot(Ax, NaN, NaN, '--r', 'LineWidth', 0.75);
% H.span = plot(Ax, [NaN, NaN], [2, 2], '-k', 'LineWidth', 2);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
title('x', 'FontSize', 10);
xlabel('Frequency (Hz)', 'FontSize', 10);
ylabel('Power (a.u.)', 'FontSize', 10);

GRP = struct();
clc
for i = 1:length(Files)
    kv = filename2struct(Files(i).name);
    idx_phen = strcmpi(PHEN.folder_name, kv.sub);
    cond = lower(PHEN.condition{idx_phen});
    EEG = load(fullfile(Files(i).folder, Files(i).name));
    chanfile = dir(sprintf('derivatives/EEG-preproc/sub-%s/ses-*/sub-*_desc-sigmanrembout_electrodes.tsv', kv.sub));
    chanlocs = readtable(fullfile(chanfile(1).folder, chanfile(1).name), 'FileType', 'text');
    if size(chanlocs, 1) ~= 179
        sprintf('%i %s', size(chanlocs, 1), Files(i).name)
    end
    chanlocs = chanlocs(1:178, :);
    FIT = EEG.FIT;
    XData = repmat([FIT(1).freq, nan], 1, length(FIT)-1);
    YData = [FIT(1:end-1).pow; nan(1, length(FIT)-1)];
    YData = YData(:);
    if any(YData > 10)
        idx_interp = find(arrayfun(@(f) any(f.pow > 10), FIT));
        for c = 1:length(idx_interp)
            d = sqrt(...
                (chanlocs.x - chanlocs.x(idx_interp(c))).^2 + ...
                (chanlocs.y - chanlocs.y(idx_interp(c))).^2 + ...
                (chanlocs.z - chanlocs.z(idx_interp(c))).^2);
            d = 1./(d.^2);
            d(idx_interp) = 0;
            d = d ./ sum(d);

            interp_data = mean([FIT(1:end-1).mu] .* d', 2);
            FIT(idx_interp(c)).mu = interp_data;

            interp_data = mean([FIT(1:end-1).peak] .* d', 2);
            FIT(idx_interp(c)).peak = interp_data;

            d = repmat(d, 1, length(FIT(1).freq));
            interp_data = mean([FIT(1:end-1).pow] .* d', 2);
            FIT(idx_interp(c)).pow = interp_data;
        end
        YData = [FIT(1:end-1).pow; nan(1, length(FIT)-1)];
        YData = YData(:);
    end
    set(H.data, 'XData', XData, 'YData', YData)
    set(H.avdata, 'XData', FIT(end).freq, 'YData', FIT(end).pow)
    set(H.peak, 'XData', [FIT(1:end-1).mu], 'YData', [FIT(1:end-1).peak])
    Ax.Title.String = sprintf('%s %s', kv.sub, cond);
    exportgraphics(Fig, sprintf('figures/subject-level/ismfit/sub-%s_ses-%s_ismfit.png', kv.sub, cond), 'Resolution', 600);
    % Store for group averages
    if isfield(GRP, cond)
        GRP.(cond)(end+1).fit = FIT;
    else
        GRP.(cond).fit = FIT;
    end
end
%%

freq = 0:0.0033:0.1;

close all

Fig = figure('Color', 'w', 'Position', [560, 640, 290 200]);
Ax = axes('NextPlot', 'add', 'FontSize', 10, 'XTick', 0:0.01:0.1, 'YTick', 1:10, 'XGrid', 'on', 'TickLength', [0 0], 'Box', 'on', 'Color', [0.96 0.97 0.99]);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
H.data = patch(Ax, 'XData', NaN, 'YData', NaN, 'EdgeColor', 'k', 'EdgeAlpha', 0.1, 'FaceColor', 'none');
H.peak = plot(Ax, NaN, NaN, '.', 'Color', [0.5 0.5 0.5], 'MarkerFaceColor', [0.5 0.5 0.5], 'MarkerSize', 6);

title('x', 'FontSize', 10);
xlabel('Frequency (Hz)', 'FontSize', 10);
ylabel('Power (a.u.)', 'FontSize', 10);

% Average across files in each condition
YData = arrayfun(@(s) s.fit, GRP.placebo, 'UniformOutput', false);
YData = cellfun(@(c) [c(1:end-1).pow], YData, 'UniformOutput', false);
YData = cat(3, YData{:});
YData = mean(YData, 3);
YData = [YData; nan(1, size(YData, 2))];
YData = YData(:);
XData = repmat([freq, nan], 1, 178);
set(H.data, 'XData', XData, 'YData', YData)
set(H.peak, 'XData', nan, 'YData', nan);
Ax.YLim = [0, 4];
Ax.Title.String = 'placebo';
exportgraphics(Fig, 'figures/ismfit_placebo.png', 'Resolution', 600);

% Average across files in each condition
YData = arrayfun(@(s) s.fit, GRP.etc120, 'UniformOutput', false);
YData = cellfun(@(c) [c(1:end-1).pow], YData, 'UniformOutput', false);
YData = cat(3, YData{:});
YData = mean(YData, 3);
YData = [YData; nan(1, size(YData, 2))];
YData = YData(:);
XData = repmat([freq, nan], 1, 178);
set(H.data, 'XData', XData, 'YData', YData)
set(H.peak, 'XData', nan, 'YData', nan);
Ax.Title.String = 'ETC120';
Ax.YLim = [0, 4];

exportgraphics(Fig, 'figures/ismfit_etc120.png', 'Resolution', 600);



