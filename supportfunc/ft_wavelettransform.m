function POW = ft_wavelettransform(EEG, Settings, varargin)
% -------------------------------------------------------------------------
% Parse optional arguments
chansel = 'all';
for i = 1:2:length(varargin)
    switch lower(varargin{i})
        case 'chansel'
            chansel = varargin{i+1};
    end
end
% -------------------------------------------------------------------------
% Select EEG channels only
switch chansel
    case 'all'
        if ~all(strcmpi({EEG.chanlocs.type}, 'EEG'))
            EEG = pop_select(EEG, 'channel', find(strcmpi({EEG.chanlocs.type}, 'EEG')));
        end
    case '1020'
        EEG = pop_select(EEG, 'channel', {'E36', 'E21', 'E224', 'E59', 'Cz', 'E183', 'E87', 'E101', 'E153', 'E116', 'E126', 'E15'});
end
% -------------------------------------------------------------------------
% Run spectral analysis for each EEG channel in the data
POW = eeg_emptyset();
T = now;
for i = 1:EEG.nbchan
    % ---------------------------------------------------------------------
    % Display remaining time
    T = remainingTime(T, EEG.nbchan, 'simple', true);
    % ---------------------------------------------------------------------
    % Update settings for this channel and store in temp-variable
    Settings.SpectrogramSettings.ChanSel = {EEG.chanlocs(i).labels};
    tmp = EEG;
    tmp.data = tmp.data(i, :);
    tmp.chanlocs = tmp.chanlocs(i);
    tmp.nbchan = 1;
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Run Spectral analysis
    tmp = Analysis_Spectrogram(tmp, Settings, {}, {});
    % ---------------------------------------------------------------------
    % If this is the first iteration, then specify the output structure
    if i == 1
        if ~strcmpi(chansel, 'all')
            POW.setname = strrep(tmp.setname, 'desc-preproc', ['desc-sigma', chansel]);
        else
            POW.setname = strrep(tmp.setname, 'desc-preproc', 'desc-sigma');
        end
        POW.filename = [POW.setname, '.set'];
        POW.filepath = tmp.filepath;
        POW.subject = tmp.subject;
        POW.nbchan = EEG.nbchan;
        POW.trials = 1;
        POW.pnts = length(tmp.spectimes);
        POW.srate = 1/tmp.spectimes(2);
        POW.xmin = 0;
        POW.xmax = tmp.spectimes(end);
        POW.times = tmp.spectimes;
        POW.chanlocs = EEG.chanlocs;
        POW.chaninfo = EEG.chaninfo;
        POW.ref = 'averef';
        POW.event = EEG.event;
        POW.data = nan(POW.nbchan, POW.pnts);
    end
    % ---------------------------------------------------------------------
    % Store the power fluctuations (averaged across frequenciesm, and square-root transformation)
    POW.data(i, :) = sqrt(mean(tmp.specdata));
end
% -------------------------------------------------------------------------
% Adjust latencies and duration of all events
for i = 1:length(POW.event)
    POW.event(i).latency = ((POW.event(i).latency - 1) ./ EEG.srate) .* POW.srate + 1;
    POW.event(i).duration = ((POW.event(i).duration - 1) ./ EEG.srate) .* POW.srate + 1;
    if isempty(POW.event(i).id)
        POW.event(i).id = max([POW.event.id]) + 1;
    end
    if isempty(POW.event(i).is_reject)
        POW.event(i).is_reject = false;
    end
end
end
