function OUT = loadarousalepochs(Taro, min_interval, sleep_episodes, arousal_type, doFooof)

Subjects = dir('derivatives/EEG-preproc/sub-*');

OUT = struct();
OUT.HR = [];
OUT.DHR = [];
OUT.SPECDATA = [];
OUT.INDEP = table();

t = now();
for i = 1:length(Subjects)

    % Load epoched data
    File = dir(sprintf('derivatives/EEG-preproc/%s/ses-*/sub*arousal_epoch.set', Subjects(i).name));
    try
        ARO = LoadDataset(fullfile(File.folder, File.name), 'all');
    catch
        continue
    end
    
    % Check we have the same number of evetns
    if sum(contains(Taro.filename, Subjects(i).name)) ~= ARO.trials
        hasAllEvents = false;
    else
        hasAllEvents = true;
    end

    % Select events
    idx_sub = contains(Taro.filename, Subjects(i).name);
    this_t = Taro(idx_sub, :);

    if ~hasAllEvents
        error('Missing events')
%         has_ev_type = {};
%         has_ev_dur = {};
%         latency = 60.*ARO.srate + 1;
%         for j = 1:length(ARO.epoch)
%             [~, this_ev_idx] = min(abs([ARO.event.latency] - latency));
%             has_ev_type = [has_ev_type; {ARO.event(this_ev_idx).type}];
%             has_ev_dur = [has_ev_dur; {ARO.event(this_ev_idx).duration./ARO.srate}];
%             latency = latency + 75.*ARO.srate;
%         end
    end

    idx_ev = ...
        this_t.interval >= min_interval & ...
        ismember(this_t.sleepepisode, sleep_episodes) & ...
        ismember(this_t.newtype, arousal_type);

    % Select data
    times = -60:0.2:14.8;
    idx_times = [];
    for j = 1:length(times)
        [~, k] = min(abs(ARO.times - times(j)));
        idx_times = [idx_times, k]; %#ok<AGROW> 
    end
    hr = squeeze(ARO.data);
    hr_bl = mean(hr(ARO.times < -1, :), 1, 'omitnan');
    hr_bl = repmat(hr_bl, ARO.pnts, 1);
    dhr = hr - hr_bl;
    hr = hr(idx_times, idx_ev);
    dhr = dhr(idx_times, idx_ev);
    specdata = ARO.specdata(:, 2:end, :, idx_ev);
    if doFooof
        nbchan = size(specdata, 1);
        fooof = ARO.fooof(idx_ev);
        fooof = arrayfun(@(f) f.f.ap_fit, fooof, 'UniformOutput', false);
        fooof = cat(4, fooof{:});
        fooof = repmat(fooof, nbchan, 1, length(ARO.spectimes), 1);
        specdata = specdata - fooof;
    end

    % Append data
    OUT.HR = cat(2, OUT.HR, hr);
    OUT.DHR = cat(2, OUT.DHR, dhr);
    OUT.SPECDATA = cat(4, OUT.SPECDATA, specdata);
    if i == 1
        OUT.specfreqs = ARO.specfreqs;
        OUT.spectimes = ARO.spectimes;
        fname = 'derivatives/EEG-preproc/sub-r0021jb/ses-1/sub-r0021jb_ses-1_task-psg_run-1_desc-nrem_channels.tsv';
        OUT.specchans = readSidecarTSV(fname, 'channels');
        OUT.specchans(~strcmpi(OUT.specchans.type, 'EEG'), :) = [];
        OUT.chanlocs = readlocs('GSN-HydroCel-257.sfp');
        OUT.chanlocs(~ismember({OUT.chanlocs.labels}, OUT.specchans.name)) = [];
    end

    % Set dep. and indep. vars
    indep = table();
    indep.stateshift = this_t.stateshift(idx_ev);
    indep.subject = repmat({Subjects(i).name}, sum(idx_ev), 1);
    indep.type = this_t.newtype(idx_ev);
    indep.duration = this_t.duration(idx_ev);
    indep.interval = this_t.interval(idx_ev);
    indep.episode = this_t.sleepepisode(idx_ev);
    OUT.INDEP = [OUT.INDEP; indep];

    t = remainingTime(t, length(Subjects), 'simple', true);
end
end