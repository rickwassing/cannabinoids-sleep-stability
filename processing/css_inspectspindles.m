function css_inspectspindles(EEG, cfg) %#ok<*INUSD> 

try
    % ---------------------------------------------------------------------
    % Get BIDS key-values
    kv = filename2struct(EEG.setname);
    % ---------------------------------------------------------------------
    % Load spindles
    SpdFiles = dir(sprintf('derivatives/EEG-processed/sub-%s/ses-%s/sub-*boxcar.set', kv.sub, kv.ses));
    SPD = [];
    for i = 1:length(SpdFiles)
        if i == 1
            SPD = LoadDataset(fullfile(SpdFiles(i).folder, SpdFiles(i).name), 'all');
        else
            SPD(i) = LoadDataset(fullfile(SpdFiles(i).folder, SpdFiles(i).name), 'all');
        end
    end
    % ---------------------------------------------------------------------
    % Create figure
    fig = figure();
    fig.Position = [1 90 1920 890];

    colors = {'red', 'green', 'blue'};
    ai = 0;
    si = 0;

    % PLOT SPINDLES

    for i = 1:3

        si = si + 1;

        ai = ai + 1;

        clr = colors{si};

        % Create axes
        ax(ai) = axes(fig); %#ok<*AGROW,*SAGROW,*LAXES>
        ax(ai).NextPlot = 'add';

        % Plot spindle density for channel 21
        YData = zeros(1, SPD(si).pnts);
        YData(SPD(si).data(21, :) > 0.9) = 1;
        YData = [0, diff(YData)];
        YData(YData <= 0) = 0;
        YData = movsum(YData, 60*SPD(si).srate+1);
        plot(SPD(si).times / (24*60*60), YData+30, '-', 'Color', standard_colors(clr).^0.5);

        % Extract 30 second interval of high and low spindle density to compare methods
        if si == 1
            tsel = struct();
            % find the maximum spindle density
            mxdens = max(YData);
            % instances where the spindle density is quite high
            idx = find(YData >= mxdens*0.8 & YData <= mxdens);
            % take a random sample out of these instances
            thisidx = randsample(idx, 35);
            % set the time selectors
            for k = 1:35
                tsel(k).hi = round((thisidx(k) - 7.5*SPD(si).srate):(thisidx(k) + 7.5*SPD(si).srate - 1));
            end
            % instances where the spindle density is quite low
            idx = find(YData >= 3 & YData <= 5);
            % take a random sample out of these instances
            thisidx = randsample(idx, 35);
            % set the time selectors
            for k = 1:35
                tsel(k).lo = round((thisidx(k) - 7.5*SPD(si).srate):(thisidx(k) + 7.5*SPD(si).srate - 1));
            end
        end

        % Plot spindle events for channel 21
        YData = SPD(si).data(21, :);
        YData(YData < 1) = 0; % remove spindles marked by arousals
        YData = [0, diff(YData)];
        YData(YData <= 0) = nan;
        if max(YData) < 2
            YData = YData.*12 + rand(1, EEG.pnts).* 3;
        end
        plot(SPD(si).times / (24*60*60), YData+30, '.', 'MarkerSize', 3, 'Color', standard_colors(clr));

        % Plot spindle density for channel 94
        YData = zeros(1, SPD(si).pnts);
        YData(SPD(si).data(94, :) > 0.9) = 1;
        YData = [0, diff(YData)];
        YData(YData <= 0) = 0;
        YData = movsum(YData, 60*SPD(si).srate+1);
        plot(SPD(si).times / (24*60*60), YData, '-', 'Color', standard_colors(clr).^0.5);

        % Plot spindle events for channel 94
        YData = SPD(si).data(94, :);
        YData(YData < 1) = 0; % remove spindles marked by arousals
        YData = [0, diff(YData)];
        YData(YData <= 0) = nan;
        if max(YData) < 2
            YData = YData.*12 + rand(1, EEG.pnts).* 3;
        end
        plot(SPD(si).times / (24*60*60), YData, '.k', 'MarkerSize', 3, 'Color', standard_colors(clr));

        % Plot divider
        plot([SPD(si).times(1), SPD(si).times(end)] / (24*60*60), [30, 30], '-k');

        % Set axis properties
        ax(ai).Box = 'off';
        ax(ai).XGrid = 'on';
        ax(ai).YLim = [0, 60];
        ax(ai).YTick = [11 , 16, 41, 46];
        ax(ai).YTickLabel = {'11', '16', '11', '16'};
        ax(ai).YGrid = 'on';
        ax(ai).OuterPosition = [0 1-si*3/15 1 3/15];
        ax(ai).TickLength = [0, 0];

        % Create topoplot
        ai = ai + 1;
        ax(ai) = axes(fig);
        skv = filename2struct(SPD(si).setname);
        topoplot(rand(1, EEG.nbchan), EEG.chanlocs, ...
            'style', 'blank', ...
            'electrodes', 'on', ...
            'hcolor', standard_colors(clr), ...
            'whitebk', 'on', ...
            'plotchans', [21 94]);
        ax(ai).Position = [0 1-si*3/15+0.05 0.1 3/15-0.1];
        ax(ai).PlotBoxAspectRatio = [1, 1, 1];
        ax(ai).XLim = [-0.8, 0.8];
        ax(ai).YLim = [-0.8, 0.8];
        ax(ai).Title.String = skv.desc;
        ax(ai).Title.Color = standard_colors(clr);

    end

    % PLOT HYPNOGRAM

    ai = ai + 1;
    ax(ai) = axes(fig);
    ax(ai).NextPlot = 'add';
    hyp = plotHypnogram(ax(ai), EEG);
    htsel(1) = plot(EEG.times([tsel(1).lo(1), tsel(1).lo(1)]) / (24*60*60), [-1.5 1.5], '-k', 'LineWidth', 3);
    htsel(2) = plot(EEG.times([tsel(1).hi(1), tsel(1).hi(1)]) / (24*60*60), [-1.5 1.5], '-k', 'LineWidth', 3);
    ax(ai).TickLength = [0, 0];
    ax(ai).XGrid = 'on';
    ax(ai).OuterPosition = [0 4/15 1 2/15];

    % PLOT EEG DATA AT SELECTED TIME SEGMENTS

    ai = ai + 1;
    ax(ai) = axes(fig);
    ax(ai).NextPlot = 'add';
    plot(EEG.times./(24*60*60), EEG.data(21, :) + 100, '-k')
    plot(EEG.times./(24*60*60), EEG.data(94, :), '-k')
    for i = 1:3
        YData = SPD(i).data(21, :);
        YData(YData < 0.9) = nan;
        YData(YData >= 0.9) = 1;
        plot(ax(ai), EEG.times./(24*60*60), YData + 100 + 35 + i * 6, '-', 'LineWidth', 3, 'Color', standard_colors(colors{i}))
        YData = SPD(i).data(94, :);
        YData(YData < 0.9) = nan;
        YData(YData >= 0.9) = 1;
        plot(ax(ai), EEG.times./(24*60*60), YData + 35 + i * 6, '-', 'LineWidth', 3, 'Color', standard_colors(colors{i}))
    end
    ax(ai).TickLength = [0, 0];
    ax(ai).XGrid = 'on';
    ax(ai).YTick = [0, 100];
    ax(ai).YTickLabel = {'0', '0'};
    ax(ai).YLim = [-100, 200];
    ax(ai).XLim = [min(EEG.times./(24*60*60)), max(EEG.times./(24*60*60))];
    ax(ai).XTick = max(EEG.times./(24*60*60));
    ax(ai).XTickLabel = {datestr(max(EEG.times./(24*60*60)), 'HH:MM')}; %#ok<*DATST>
    ax(ai).OuterPosition = [0 0 0.5 4/15];
    ax(ai).Position(1) = 0.01;
    ax(ai).Position(3) = 0.48;

    ai = ai + 1;
    ax(ai) = axes(fig);
    ax(ai).NextPlot = 'add';
    plot(EEG.times./(24*60*60), EEG.data(21, :) + 100, '-k')
    plot(EEG.times./(24*60*60), EEG.data(94, :), '-k')
    for i = 1:3
        YData = SPD(i).data(21, :);
        YData(YData < 0.9) = nan;
        YData(YData >= 0.9) = 1;
        plot(ax(ai), EEG.times./(24*60*60), YData + 100 + 35 + i * 6, '-', 'LineWidth', 3, 'Color', standard_colors(colors{i}))
        YData = SPD(i).data(94, :);
        YData(YData < 0.9) = nan;
        YData(YData >= 0.9) = 1;
        plot(ax(ai), EEG.times./(24*60*60), YData + 35 + i * 6, '-', 'LineWidth', 3, 'Color', standard_colors(colors{i}))
    end
    ax(ai).TickLength = [0, 0];
    ax(ai).XGrid = 'on';
    ax(ai).YTick = [0, 100];
    ax(ai).YTickLabel = {'0', '0'};
    ax(ai).YLim = [-100, 200];
    ax(ai).XLim = [min(EEG.times./(24*60*60)), max(EEG.times./(24*60*60))];
    ax(ai).XTick = max(EEG.times./(24*60*60));
    ax(ai).XTickLabel = {datestr(max(EEG.times./(24*60*60)), 'HH:MM')}; %#ok<*DATST>
    ax(ai).OuterPosition = [0.5 0 0.5 4/15];
    ax(ai).Position(1) = 0.51;
    ax(ai).Position(3) = 0.48;

    % ALIGN AXES

    ax(1).XLim = ax(7).XLim;
    ax(1).XTick = ax(7).XTick;
    ax(1).XTickLabel = {};

    ax(3).XLim = ax(7).XLim;
    ax(3).XTick = ax(7).XTick;
    ax(3).XTickLabel = {};

    ax(5).XLim = ax(7).XLim;
    ax(5).XTick = ax(7).XTick;
    ax(5).XTickLabel = {};

    % Creat output directory
    outfilepath = sprintf('inspect/spindledet/sub-%s/ses-%s', kv.sub, kv.ses);
    if exist(outfilepath, 'dir') == 0
        mkdir(outfilepath)
    end

    for i = 1:length(tsel)
        ax(end-1).XLim = [min(EEG.times(tsel(i).lo)./(24*60*60)), max(EEG.times(tsel(i).lo)./(24*60*60))];
        ax(end).XLim = [min(EEG.times(tsel(i).hi)./(24*60*60)), max(EEG.times(tsel(i).hi)./(24*60*60))];
        ax(end-1).XTick = max(EEG.times(tsel(i).lo)./(24*60*60));
        ax(end).XTick = max(EEG.times(tsel(i).hi)./(24*60*60));
        ax(end-1).XTickLabel = {datestr(max(EEG.times(tsel(i).lo)./(24*60*60)), 'HH:MM')}; %#ok<*DATST>
        ax(end).XTickLabel = {datestr(max(EEG.times(tsel(i).hi)./(24*60*60)), 'HH:MM')}; %#ok<*DATST>
        htsel(1).XData = EEG.times([tsel(i).lo(1), tsel(i).lo(1)]) / (24*60*60);
        htsel(2).XData = EEG.times([tsel(i).hi(1), tsel(i).hi(1)]) / (24*60*60);

        fullfilename = sprintf('%s/sub-%s_ses-%s_spindle-%i_inspect.png', outfilepath, kv.sub, kv.ses, i);
        exportgraphics(fig, fullfilename, 'Resolution', 144);
    end

    close all

catch ME %#ok<*NASGU> 
    keyboard
end

end
