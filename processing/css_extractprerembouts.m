function css_extractprerembouts(filepath, cfg)
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
% Extract hypnogram information
Thyp = css_eeglab2hypnogram(EEG);
dur = [asrow(find(diff([nan; Thyp.episode == 2]) == 1)); asrow(find(diff([nan; Thyp.episode == 2]) == -1))];
dur = diff(dur) * 30 * EEG.srate; % convert to samples
lat = ascolumn(find(diff([nan; Thyp.episode == 2]) == 1)-1) * 30 * EEG.srate + 1; % convert to samples
bouts = ascolumn(find(diff([nan; Thyp.episode == 2]) == 1));
bouts = bouts * 30 - 30; % Convert to seconds
bouts = [bouts-300, bouts+60];
% Insert episode events into the EEG.event struct
for i = 1:size(bouts, 1)
    EEG.event(end+1).latency = lat(i);
    EEG.event(end).duration = dur(i);
    EEG.event(end).type = 'remeps';
    EEG.event(end).id = max([EEG.event.id])+1;
    EEG.event(end).is_reject = false;
end
% -------------------------------------------------------------------------
% Save original event latencies
EEG.event = storeoriglatency(EEG.event);
EEG = eeg_checkset(EEG, 'eventconsistency');
% -------------------------------------------------------------------------
% Select EEG during N2 and N3
BOUT = pop_select(EEG, 'time', bouts);
BOUT = eeg_checkset(BOUT, 'eventconsistency');
BOUT.times = BOUT.xmin:1/BOUT.srate:BOUT.xmax; % undo 'pop_select' convertion to milliseconds
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Save
[BOUT.filepath, BOUT.setname] = fileparts(cfg.outfilepath);
BOUT.filename = [BOUT.setname, '.set'];
BOUT = SaveDataset(BOUT, 'all');
% -------------------------------------------------------------------------
Fig = plotnremboutselection(EEG, bouts);
exportgraphics(Fig, fullfile(BOUT.filepath, [BOUT.setname, '.png']), 'Resolution', 600)

end
