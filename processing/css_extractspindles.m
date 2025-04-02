function css_extractspindles(filepath, cfg)
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
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Only keep EEG channels
    if ~all(strcmpi({EEG.chanlocs.type}, 'EEG'))
        EEG = pop_select(EEG, 'channel', find(strcmpi({EEG.chanlocs.type}, 'EEG')));
    end
    % -------------------------------------------------------------------------
    % Get hypnogram table
    HYP = eeglab2hypnogram(EEG);
    try
        HYP.sleepstage = HYP.stage_num;
    catch ME
        keyboard
    end
    % Generate vector with as many data points as the EEG data and store the Hypnogram
    v_Hyp = zeros(1, EEG.pnts, 'single');
    for i = 1:size(HYP, 1)
        idx = round(HYP.latency(i)+1);
        idx = idx:(idx+30*EEG.srate-1);
        switch HYP.stage_str{i}
            case 'wake'
                v_Hyp(idx) = 0;
            case 'n1'
                v_Hyp(idx) = -1;
            case 'n2'
                v_Hyp(idx) = -2;
            case 'n3'
                v_Hyp(idx) = -3;
            case 'rem'
                v_Hyp(idx) = -4;
        end
    end
    v_Hyp = v_Hyp(1:EEG.pnts);
    % -------------------------------------------------------------------------
    % Get vector of arousal events to exclude from spindle detection
    st_AroEvents = EEG.event(ismember({EEG.event.type}, {'arousal', 'arousalemg'}));
    v_Aro = events2timeseries(st_AroEvents, EEG.xmin, EEG.xmax, EEG.srate);
    % -------------------------------------------------------------------------
    % Create timeseries of spindle activity
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Detect spindles and store the cumulative sum in 5 second moving windows
    rt = now();
    for i = 1:EEG.nbchan
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        % Detect spindles
        switch lower(cfg.method)
            case 'fernandez'
                [m_Spindles, st_Spindles] = f_SpDetection_Humans(...
                    double(EEG.data(i, :)), ...
                    EEG.srate, ...
                    v_Hyp, EEG, i);
                v_SpindleFreqs = ones(1, size(m_Spindles, 2));
            case 'ferrarelli'
                [m_Spindles, ~, v_SpindleFreqs] = f_SpDetection_Ferrarelli(...
                    double(EEG.data(i, :)), ...
                    EEG.srate, ...
                    v_Hyp, EEG, i);
            case 'wamsley'
                [m_Spindles, ~, v_SpindleFreqs] = f_SpDetection_Wamsley(...
                    double(EEG.data(i, :)), ...
                    EEG.srate, ...
                    v_Hyp, EEG, i);
        end
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        % Create spindle timeseries signal
        v_SpSignal = zeros(1, EEG.pnts);
        for j = 1:size(m_Spindles, 2)
            v_SpSignal(m_Spindles(1, j):m_Spindles(2, j)) = v_SpindleFreqs(j);
        end
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        % The v_SpSignal contains non-zero values for each spindle event,
        % where the value indicates the mean frequency of the spindle. Now
        % mutiply these values by -1 for each arousal event. This way we
        % can inspect instances where an arousal was scored and a spindle
        % was detected (false positive spindle, or false positive arousal)
        v_Mult = ones(size(v_SpSignal));
        v_Mult(v_Aro) = -1;
        v_SpSignal = v_SpSignal .* v_Mult;
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        % Store in data structure
        EEG.data(i, :) = single(v_SpSignal);
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        % Show remaining time
        rt = remainingTime(rt, EEG.nbchan, true);
    end
    % ---------------------------------------------------------------------
    % Save dataset
    [EEG.filepath, EEG.setname] = fileparts(cfg.outfilepath);
    EEG.filename = [EEG.setname, '.set'];
    SaveDataset(EEG, 'all');
    
catch ME
    keyboard
end
end
