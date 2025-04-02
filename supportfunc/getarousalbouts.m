function [EEG, arousalbouts] = getarousalbouts(EEG, HYP, cfg)

% Legacy support
if ~any(strcmpi(HYP.Properties.VariableNames, 'sleepstage'))
    HYP.sleepstage = HYP.stage_num;
end
arousalbouts = [];
idx_sel = [];
for i = 1:length(EEG.event)
    % ---------------------------------------------------------------------
    % Init the event 'stage' field
    EEG.event(i).stage = NaN;
    EEG.event(i).next_stage = NaN;
    EEG.event(i).is_selected = false;
    EEG.event(i).is_awakening = false;
    % ---------------------------------------------------------------------
    % If it is not an arousal event, then continue
    if ~contains(EEG.event(i).type, {'alpha', 'arousal'})
        continue
    end
    % ---------------------------------------------------------------------
    % Init
    is_req_stage = false; % boolean to check if the arousal occurred in the requested sleep stage
    % ---------------------------------------------------------------------
    % Find the index of the epoch this arousal occurred in
    latency = EEG.event(i).latency;
    idx_epoch = find(HYP.latency <= latency, 1, 'last');
    epoch_latency = (latency - HYP.latency(idx_epoch)) / EEG.srate;
    % ---------------------------------------------------------------------
    % Check the current sleep stage and previous one
    this_stage = HYP.sleepstage(idx_epoch);
    if idx_epoch > 1
        prev_stage = HYP.sleepstage(idx_epoch-1);
    else
        prev_stage = 1;
    end
    if idx_epoch < length(HYP.sleepstage)
        next_stage = HYP.sleepstage(idx_epoch+1);
    else
        next_stage = 1;
    end
    % ---------------------------------------------------------------------
    % Check if this is an arousal associated with the requested 'stage'
    if this_stage == cfg.stage
        % Arousal onset occured in the requested stage
        is_req_stage = true;
        EEG.event(i).stage = this_stage;
    elseif epoch_latency <= 15 && prev_stage == cfg.stage
        % Arousal onset occurred in the first 15 seconds of an epoch that
        % was preceeded by the requested stage
        is_req_stage = true;
        EEG.event(i).stage = prev_stage;
    end
    % ---------------------------------------------------------------------
    % Store the stage as an event field
    EEG.event(i).next_stage = next_stage;
    % Store if the arousal resulted in an awakening or not
    if EEG.event(i).stage <= -2 || EEG.event(i).stage == 0
        if next_stage == -1 || next_stage == 1
            EEG.event(i).is_awakening = true;
        end
    end
    % ---------------------------------------------------------------------
    % Continue if it is not an associated arousal
    if ~is_req_stage
        continue
    end
    % ---------------------------------------------------------------------
    % Check if at least 'cutoff' seconds of pre-arousal sleep was the requested 'stage'
    cnt = 0;
    bout_start = idx_epoch - 1;
    while bout_start > 0
        % The current duration of the selected bout
        bout_dur = latency - HYP.latency(bout_start);
        % If this bout-start epoch is N2, reset the counter
        if HYP.sleepstage(bout_start) == cfg.stage
            cnt = 0;
        end
        % If this bout-start epoch is not N2, increment the counter
        if HYP.sleepstage(bout_start) ~= cfg.stage
            cnt = cnt+1;
        end
        % If there are more than one consecutive not-'stage' epochs, then break
        if cnt > 1
            bout_start = bout_start + 2; % revert last seen epoch that was 'stage'
            bout_dur = latency - HYP.latency(bout_start);
            break
        end
        % If we checked at least 'cutoff' seconds, we can stop
        if bout_dur >= cfg.cutoff*EEG.srate
            break
        end
        % In the next iteration, check the epoch before the current one
        bout_start = bout_start - 1;
    end
    % ---------------------------------------------------------------------
    % Check the length of the pre-arousal bout
    if bout_dur >= (cfg.cutoff*EEG.srate)
        % Add this bout to the selection. The onset of the bout must be an
        % integer number of seconds pre-arousal and the offset of the bout 
        % is 30 seconds after the onset of the arousal so we can check if 
        % the arousal was followed by an awakening or continued sleep.
        whole_seconds = latency:-EEG.srate:HYP.latency(bout_start);
        this_bout_start = whole_seconds(end);
        arousalbouts = [arousalbouts; this_bout_start, round(latency+30*EEG.srate)-1]; %#ok<AGROW> 
        % Store the index of the selected events
        idx_sel = [idx_sel; i]; %#ok<AGROW> 
    end
end
% Store which events have been selected in the event struct
for i = 1:length(idx_sel)
    EEG.event(idx_sel(i)).is_selected = true;
end
end
