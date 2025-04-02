function css_extractnrembouts(filepath, cfg)
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
% Save original event latencies
EEG.event = storeoriglatency(EEG.event);
% -------------------------------------------------------------------------
% Extract hypnogram information
Thyp = css_eeglab2hypnogram(EEG);
% -------------------------------------------------------------------------
% Select EEG during N2 and N3
nrembouts = getnrembouts(Thyp, EEG.srate, 300); % '300' is minimum bout duration in seconds
NREM = pop_select(EEG, 'time', nrembouts);
NREM = eeg_checkset(NREM, 'eventconsistency');
NREM.times = NREM.xmin:1/NREM.srate:NREM.xmax; % undo 'pop_select' convertion to milliseconds
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Save
[NREM.filepath, NREM.setname] = fileparts(cfg.outfilepath);
NREM.filename = [NREM.setname, '.set'];
NREM = SaveDataset(NREM, 'all');
% -------------------------------------------------------------------------
Fig = plotnremboutselection(EEG, nrembouts);
exportgraphics(Fig, fullfile(NREM.filepath, [NREM.setname, '.png']), 'Resolution', 600)

end
