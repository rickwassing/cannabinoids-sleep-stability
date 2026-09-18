function ARO = css_extractarousalbouts(filepath, cfg)
% -------------------------------------------------------------------------
% Default values
def.allowalpha = true;
def.allowmicro = true;
fnames = fieldnames(cfg);
for i = 1:length(fnames)
    def.(fnames{i}) = cfg.(fnames{i});
end
cfg = def;
% check config
reqfields = {'stage', 'cutoff'};
for i = 1:length(reqfields)
    if ~isfield(cfg, reqfields{i})
        error('Configuration must contain the field ''%s''.', reqfields{i})
    end
end
% -------------------------------------------------------------------------
% Check if the 'filepath' variable is a string (path to file) or a stuct
% (preloaded data in EEGLAB structure)
if isstruct(filepath)
    EEG = filepath; % 'EEG' struct was used instead of path
    clear filepath;
else
    EEG = css_loadeegandfilter(filepath);
end
% -------------------------------------------------------------------------
% Keep only some events
if cfg.allowalpha
    allowedevents = {'wake', 'n1', 'n2', 'n3', 'rem', 'alpha', 'alphaemg', 'arousal', 'arousalemg', 'loff', 'lon', 'reject'};
else
    allowedevents = {'wake', 'n1', 'n2', 'n3', 'rem', 'arousal', 'arousalemg', 'loff', 'lon', 'reject'};
end
idx_keep = ismember({EEG.event.type}, allowedevents);
EEG.event = EEG.event(idx_keep);
% Round the onset latency to the nearest integer
for i = 1:length(EEG.event)
    EEG.event(i).latency = round(EEG.event(i).latency);
end
% -------------------------------------------------------------------------
% Save original event latencies
EEG.event = storeoriglatency(EEG.event);
% -------------------------------------------------------------------------
% Rename arousals shorter than 3 seconds to microarousals
idx_rm = [];
for i = 1:length(EEG.event)
    if ~ismember(EEG.event(i).type, {'alpha', 'alphaemg', 'arousal', 'arousalemg'})
        continue
    end
    if EEG.event(i).duration < 3*EEG.srate
        EEG.event(i).type = ['micro', EEG.event(i).type];
        idx_rm = [idx_rm; i]; %#ok<AGROW>
    end
end
% remove micro-arousals if they are not allowed
if ~cfg.allowmicro
    EEG.event(idx_rm) = [];
end
% -------------------------------------------------------------------------
% Select arousal events that occurred in association with the requested
% sleep stage. That is, the arousal onset occurred in the requested stage
% or within the first 15 seconds of an epoch following the requested stage.
% -------------------------------------------------------------------------
% Get hypnogram table
HYP = css_eeglab2hypnogram(EEG);
% -------------------------------------------------------------------------
% Find the arousal bouts
[EEG, arousalbouts] = getarousalbouts(EEG, HYP, cfg);
keyboard
% -------------------------------------------------------------------------
% Select parts of the data for each arousal and their pre-arousal N2 sleep
if ~isempty(arousalbouts)
    for i = 1:size(arousalbouts, 1)
        fprintf('Extracing bout %i of %i (%.0f%%)\n', i, size(arousalbouts, 1), 100*i/size(arousalbouts, 1))
        if i == 1
            ARO = pop_select(EEG, 'point', arousalbouts(i, :));
            ARO.event(end+1).latency = ARO.pnts-0.5;
            ARO.event(end).duration = 0;
            ARO.event(end).type = 'boundary';
            ARO.event(end).id = max([EEG.event.id])+1;
            ARO.event(end).is_reject = false;
            ARO.event(end).origlatency = ARO.pnts-0.5;
            ARO.event(end).stage = nan;
            ARO.event(end).next_stage = nan;
            ARO.event(end).is_selected = false;
            ARO.event(end).is_awakening = false;
        else
            tmp = pop_select(EEG, 'point', arousalbouts(i, :));
            tmp.event(end+1).latency = tmp.pnts-0.5;
            tmp.event(end).duration = 0;
            tmp.event(end).type = 'boundary';
            tmp.event(end).id = max([ARO.event(end).id])+1;
            tmp.event(end).is_reject = false;
            tmp.event(end).origlatency = tmp.pnts-0.5;
            tmp.event(end).stage = nan;
            tmp.event(end).next_stage = nan;
            tmp.event(end).is_selected = false;
            tmp.event(end).is_awakening = false;
            for j = 1:length(tmp.event)
                tmp.event(j).latency = tmp.event(j).latency + ARO.pnts;
            end            
            ARO.data = [ARO.data, tmp.data];
            ARO.event = [ARO.event, tmp.event];
            ARO.pnts = size(ARO.data, 2);
            ARO.xmax = (ARO.pnts-1)/ARO.srate;
            ARO.times = 0:1/ARO.srate:ARO.xmax;
        end
    end
    ARO = eeg_checkset(ARO, 'eventconsistency');
    ARO = forceValidEventType(ARO);
    ARO.times = ARO.xmin:1/ARO.srate:ARO.xmax; % undo 'pop_select' convertion to milliseconds
else
    return; % No arousal bouts were selected...
end
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Save
if nargout == 0
    [ARO.filepath, ARO.setname] = fileparts(cfg.outfilepath);
    ARO.filename = [ARO.setname, '.set'];
    ARO = SaveDataset(ARO, 'all');
    % -------------------------------------------------------------------------
    close all;
    Fig = plotarousalboutselection(EEG);
    exportgraphics(Fig, fullfile(ARO.filepath, [ARO.setname, '.png']), 'Resolution', 600)
end
end
