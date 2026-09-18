function ARO = eeglab2arousals(EEG, HYP)
ARO = EEG.event(ismember(lower({EEG.event.type}), {'alpha', 'alphaemg', 'arousal', 'arousalemg'}));
if isempty(ARO)
    return
end
for i = 1:length(ARO)
    idx = find(HYP.latency <= ARO(i).latency, 1, 'last');
    if isempty(idx)
        continue
    end
    if i == 1
        ARO(i).interval = nan;
    else
        ARO(i).interval = ARO(i).latency - (ARO(i-1).latency + ARO(i-1).duration);
    end
    ARO(i).currsleepstage = HYP.sleepstage(idx);
    ARO(i).sleepcycle = HYP.sleepcycle(idx);
    ARO(i).sleepepisode = HYP.sleepepisode(idx);
    aroinepoch = [ARO([ARO.latency] >= HYP.latency(idx) & [ARO.latency] <  HYP.latency(idx)+30*EEG.srate).id];
    if isempty(aroinepoch)
        ARO(i).isfirstinepoch = false;
    elseif ARO(i).id == aroinepoch(1)
        ARO(i).isfirstinepoch = true;
    else
        ARO(i).isfirstinepoch = false;
    end
    if isempty(aroinepoch)
        ARO(i).islastinepoch = false;
    elseif ARO(i).id == aroinepoch(end)
        ARO(i).islastinepoch = true;
    else
        ARO(i).islastinepoch = false;
    end
    if  idx > 1
        ARO(i).prevsleepstage = HYP.sleepstage(idx-1);
    else
        ARO(i).prevsleepstage = nan;
    end
    if idx < length(HYP.sleepstage)
        ARO(i).nextsleepstage = HYP.sleepstage(idx+1);
    else
        ARO(i).nextsleepstage = nan;
    end
end
ARO = struct2table(ARO);
try
    ARO = movevars(ARO, 'interval', 'After', 'duration');
catch ME
    keyboard
end
ARO = movevars(ARO, 'prevsleepstage', 'Before', 'currsleepstage');
ARO = movevars(ARO, 'nextsleepstage', 'After', 'currsleepstage');
ARO.epoch = ceil((ARO.latency / EEG.srate)/30);
% If this arousal occurred in W or in N1 where the previous epoch was sleep, and it is the first arousal, and the previous epoch did not have a state-shift arousal, then it is a state-shift arousal
% If this arousal occurred in W or in N1 where the previous epoch was sleep, and it is the first arousal, and the previous epoch did have a state-shift arousal, then it is not a state-shift arousal
% If this arousal occurred in sleep (N2, N3, R), and it is the last arousal in the epoch, and the next epoch is W or N1, then it is a state-shift arousal
% If this arousal occurred in sleep (N2, N3, R), and it is not the last arousal, it is not a state-shift arousal
ARO.stateshift = false(size(ARO, 1), 1);
clamp = false; % Assume a state shift can occur
clampoffset = nan;
clampfrom = nan;
for i = 1:size(ARO, 1)
    % Check if we can release the clamp i.e., if at least one epoch of
    % sleep occurred since the last state shift
    if clamp
        idx = find(HYP.latency <= clampoffset, 1, 'last'):find(HYP.latency <= ARO.latency(i), 1, 'last');
        if clampfrom == -1
            if any(ismember(HYP.sleepstage(idx), -1))
                clamp = false;
            end
        else
            if ...
                    any(ismember(HYP.sleepstage(idx), [0, -2, -3])) || ...
                    ~isempty(strfind(asrow(HYP.sleepstage(idx)), [-1, -1])) || ...
                    ~isempty(strfind(asrow(HYP.sleepstage(idx)), [1, -1])) %#ok<STREMP>
                clamp = false;
            end
        end

    end
    % If clamp is true, then we don't have to check for state shifts
    if clamp
        ARO.stateshift(i) = false;
        ARO.shiftfrom(i) = nan;
        ARO.shiftto(i) = nan;
        continue
    end
    % Check if a state shift occurred
    if i == 1
        ARO.stateshift(i) = ...
            (ARO.currsleepstage(i) == 1 & ARO.prevsleepstage(i) < 1 & ARO.isfirstinepoch(i)) | ...
            (ARO.currsleepstage(i) == -1 & ARO.prevsleepstage(i) < 1 & ARO.prevsleepstage(i) ~= -1 & ARO.isfirstinepoch(i)) | ...
            (ARO.currsleepstage(i) == 0 | ARO.currsleepstage(i) <= -2) & ARO.islastinepoch(i) & (ARO.nextsleepstage(i) == -1 | ARO.nextsleepstage(i) == 1);
    else
        idx = ARO.epoch == ARO.epoch(i)-1;
        if any(idx)
            prevContainsShift = any(ARO.stateshift(idx));
        else
            prevContainsShift = false;
        end
        shift_inwake = (ARO.currsleepstage(i) == 1 & ARO.prevsleepstage(i) < 1 & ARO.isfirstinepoch(i) & ~prevContainsShift);
        shift_inn1 = (ARO.currsleepstage(i) == -1 & ARO.prevsleepstage(i) < 1 & ARO.prevsleepstage(i) ~= -1 & ARO.isfirstinepoch(i) & ~prevContainsShift);
        shift_fromn1 = (ARO.currsleepstage(i) == -1 & ARO.islastinepoch(i) & ARO.nextsleepstage(i) == 1 & ~prevContainsShift);
        shift_fromsleep = (ARO.currsleepstage(i) == 0 | ARO.currsleepstage(i) <= -2) & ARO.islastinepoch(i) & (ARO.nextsleepstage(i) == -1 | ARO.nextsleepstage(i) == 1);
        ARO.stateshift(i) = shift_inwake | shift_inn1 | shift_fromn1 | shift_fromsleep;
    end
    % If a state shift occurred, save the from and to
    if ~ARO.stateshift(i)
        ARO.shiftfrom(i) = nan;
        ARO.shiftto(i) = nan;
    else
        if shift_inwake || shift_inn1
            ARO.shiftfrom(i) = ARO.prevsleepstage(i);
            ARO.shiftto(i) = ARO.currsleepstage(i);
        elseif shift_fromn1 || shift_fromsleep
            ARO.shiftfrom(i) = ARO.currsleepstage(i);
            ARO.shiftto(i) = ARO.nextsleepstage(i);
        end
    end
    % If a state shift occurred, set clamp to true, until one epoch of sleep has occurred
    if ARO.stateshift(i)
        clamp = true;
        clampoffset = ARO.latency(i) + ARO.duration(i);
        clampfrom = ARO.shiftfrom(i);
    end
end

end