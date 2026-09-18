function css_detectarousals(filepath, cfg)
try
    % -------------------------------------------------------------------------
    % Check if the 'filepath' variable is a string (path to file) or a stuct
    % (preloaded data in EEGLAB structure)
    if isstruct(filepath)
        EEG = filepath; % 'EEG' struct was used instead of path
        clear filepath;
    else
        EEG = css_loadeegandfilter(filepath);
    end

    % We also need EMG, which is not in the preproc file (silly I excluded
    % all the phys channels when saving the preproc files, but that's what
    % I have now, so will just load this here
    PHEN = readtable('Scoring log and notes_CANSLEEP arousal.xlsx');
    kv = filename2struct(EEG.setname);
    idx_phen = find(strcmpi(PHEN.r_number, kv.sub) & strcmpi(PHEN.condition, kv.ses));
    if isempty(idx_phen) || length(idx_phen) > 1
        keyboard
    end
    filepath = dir(sprintf('./derivatives/EEG-inspect/sub-%s/sub-*desc-inspect_eeg.set', PHEN.folder_name{idx_phen}));
    if isempty(filepath) || length(filepath) > 1
        keyboard
    end
    PHYS = LoadDataset(fullfile(filepath.folder, filepath.name), 'all');
    idx_emg = find(strcmpi({PHYS.chanlocs.type}, 'emg'));
    PHYS = pop_select(PHYS, 'channel', idx_emg);
    % Apply filter
    fcutoff = 15;
    forder = pop_firwsord('hamming', PHYS.srate, 1);
    PHYS = pop_firws(PHYS, ...
        'fcutoff', fcutoff, ...
        'ftype', 'highpass', ...
        'wtype', 'hamming', ...
        'forder', forder, ...
        'minphase', 0);
    for i = 1:size(PHYS.data, 1)
        PHYS.data(i, PHYS.data(i, :) > 2000) = 2000;
        PHYS.data(i, PHYS.data(i, :) < -2000) = -2000;
    end
    PHYS.data = abs(PHYS.data);
    % 200 ms moving average
    PHYS.data = movmean(PHYS.data, 0.2*PHYS.srate+1, 2);
    % Resample
    PHYS = pop_resample(PHYS, 128);

    % Get channel indices
    idx_fz = find(strcmpi({EEG.chanlocs.labels}, 'E21')); %#ok<NASGU> 
    idx_cz = find(strcmpi({EEG.chanlocs.labels}, 'Cz'));
    idx_pz = find(strcmpi({EEG.chanlocs.labels}, 'E101')); %#ok<NASGU> 

    % Construct the EEG and EMG signal
    eeg = struct();
    eeg.raw = ascolumn(EEG.data(idx_cz, :));
    eeg.rate = EEG.srate;
    
    emg = struct();
    emg.raw = ascolumn(mean(PHYS.data(1:2, :)));
    emg.rate = PHYS.srate;

    % Run detect arousals
    [arousals, ~, ~, ~] = arousalDetection(eeg, emg);

    % Remove arousals that occurred prior to fist sleep epoch and after last sleep epoch
    idx_son = find(ismember({EEG.event.type}, {'n1', 'n2', 'n3', 'rem'}), 1, 'first');
    idx_soff = find(ismember({EEG.event.type}, {'n1', 'n2', 'n3', 'rem'}), 1, 'last');
    idx_rm = [arousals.endSample] < EEG.event(idx_son).latency;
    arousals(idx_rm) = [];
    idx_rm = [arousals.startSample] > EEG.event(idx_soff).latency;
    arousals(idx_rm) = [];

    % Remove arousals that were not preceeded by 10 seconds of sleep
    idx_rm = [];
    for i = 1:length(arousals)
        % did this arousal occur in a sleep epoch, or a wake epoch
        % preceeded by a sleep epoch
        idx_check = find([EEG.event.latency] < arousals(i).startSample & ismember({EEG.event.type}, {'n1', 'n2', 'n3', 'rem', 'wake'}), 2, 'last');
        if all(strcmpi({EEG.event(idx_check).type}, 'wake'))
            idx_rm = [idx_rm, i]; %#ok<AGROW> 
        end
    end
    arousals(idx_rm) = [];

    % Add detected arousals 
    for i = 1:length(arousals)
        label = 'auto_';
        switch arousals(i).source
            case 'alpha'
                label = [label, 'alpha']; %#ok<AGROW> 
            otherwise
                label = [label, 'arousal']; %#ok<AGROW> 
        end
        if any(strcmpi(arousals(i).info, 'emg activity'))
            label = [label, 'emg']; %#ok<AGROW> 
        end
        EEG.event(end+1).latency = arousals(i).startSample;
        EEG.event(end).duration = arousals(i).duration*128;
        EEG.event(end).type = label;
        EEG.event(end).id = max([EEG.event.id])+1;
        EEG.event(end).is_reject = false;
    end
    EEG = eeg_checkset(EEG, 'eventconsistency');

    idx_stage = find(ismember({EEG.event.type}, {'n1', 'n2', 'n3', 'rem', 'wake'}));
    idx_aro = find(contains({EEG.event.type}, 'alpha') | contains({EEG.event.type}, 'arousal'));

    close all;
    Fig = figure('Color', 'w', 'Position', [1 600 1920 375]);
    Ax(1) = axes(Fig, ...
        'Box', 'on', ...
        'TickLength', [0 0], ...
        'OuterPosition', [0 0.5 1 0.5], ...
        'NextPlot', 'add');
    plot(Ax(1), EEG.times, EEG.data(idx_cz, :), '-k')
    Ax(2) = axes(Fig, ...
        'Box', 'on', ...
        'TickLength', [0 0], ...
        'OuterPosition', [0 0 1 0.5], ...
        'NextPlot', 'add');
    plot(Ax(2), EEG.times, emg.raw, '-k')

    for i = 1:length(idx_stage)
        XData = [...
            EEG.event(idx_stage(i)).latency; ...
            EEG.event(idx_stage(i)).latency];
        YData = [-150 150];
        plot(Ax(1), XData./EEG.srate, YData, ':k')
        plot(Ax(2), XData./EEG.srate, YData, ':k')
        text(Ax(1), XData(1)./EEG.srate, YData(1), sprintf(' %s', EEG.event(idx_stage(i)).type), 'FontSize', 10, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom')
    end

    for i = 1:length(idx_aro)
        XData = [...
            EEG.event(idx_aro(i)).latency; ...
            EEG.event(idx_aro(i)).latency + EEG.event(idx_aro(i)).duration; ...
            EEG.event(idx_aro(i)).latency + EEG.event(idx_aro(i)).duration; ...
            EEG.event(idx_aro(i)).latency; ...
            EEG.event(idx_aro(i)).latency];
        if contains(EEG.event(idx_aro(i)).type, 'auto')
            YData = [-140 -140 140 140 -140];
        else
            YData = [-100 -100 100 100 -100];
        end
        plot(Ax(1), XData./EEG.srate, YData, '--r')
        plot(Ax(2), XData./EEG.srate, YData, '--r')
        text(Ax(1), XData(1)./EEG.srate, YData(1), sprintf(' %s %i', strrep(EEG.event(idx_aro(i)).type, '_', ' '), EEG.event(idx_aro(i)).id), 'FontSize', 10, 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom')
    end
    Ax(1).YLim = [-150 150];
    Ax(2).YLim = [-100 100];

    for i = 1:length(idx_aro)
        Ax(1).XLim = [EEG.event(idx_aro(i)).latency/EEG.srate - 10, EEG.event(idx_aro(i)).latency/EEG.srate + 20];
        Ax(2).XLim = [EEG.event(idx_aro(i)).latency/EEG.srate - 10, EEG.event(idx_aro(i)).latency/EEG.srate + 20];
        filename = strrep(EEG.setname, '_eeg', sprintf('_%i%s%i.png', round(EEG.event(idx_aro(i)).latency), EEG.event(idx_aro(i)).type, EEG.event(idx_aro(i)).id));
        exportgraphics(Fig, ['inspect/arousals/', filename], 'Resolution', 150)
    end

    % ---------------------------------------------------------------------
    % Save dataset
    [EEG.filepath, EEG.setname] = fileparts(cfg.outfilepath);
    EEG.filename = [EEG.setname, '.set'];
    SaveDataset(EEG, 'header');
    
catch ME %#ok<NASGU> 
    keyboard
end
end
