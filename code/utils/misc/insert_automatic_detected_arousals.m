Files = dir('derivatives/EEG-preproc/sub-*/ses-*/sub-*desc-preproc*_eeg.set');
Evs = dir('derivatives/arousal_scoring/outputs/finetuned-batchnorm/events/sub-*.tsv');
T = readtable('Scoring log and notes_CANSLEEP arousal.xlsx');
clc
for i = 1:length(Files)
    kv = filename2struct(Files(i).name);
    if startsWith(kv.sub, 'r')
        visit = T.folder_name{strcmpi(T.r_number, kv.sub) & strcmpi(T.condition, kv.ses)};
        idx_ev = find(contains({Evs.name}, visit));
    else
        idx_ev = find(contains({Evs.name}, kv.sub));
    end

    try
        hdr = LoadDataset(fullfile(Files(i).folder, Files(i).name), 'header');
        event = readtable(fullfile(Evs(idx_ev).folder, Evs(idx_ev).name), 'FileType', 'text');
    catch
        continue
    end
    idx_rm = strcmpi({hdr.event.type}, 'detectarousal');
    hdr.event(idx_rm) = [];
    for ev = 1:size(event, 1)
        hdr.event(end+1).latency = event.onset(ev)*hdr.srate;
        hdr.event(end).duration = event.duration(ev)*hdr.srate;
        hdr.event(end).type = 'detectarousal';
        hdr.event(end).is_reject = false;
    end
    for ev = 1:length(hdr.event)
        hdr.event(ev).id = ev;
    end
    hdr = eeg_checkset(hdr, 'eventconsistency');
    if isfield(hdr.etc, 'viewsettings')
        hdr.etc = rmfield(hdr.etc, 'viewsettings');
    end
    if isfield(hdr.etc, 'faster')
        hdr.etc = rmfield(hdr.etc, 'faster');
    end
    SaveDataset(hdr, 'header')
end




%%
SaveDataset(HDR, 'header')