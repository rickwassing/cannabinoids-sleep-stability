function EEG = css_createfstlvloutput(outfname, Features, DataStr)

if nargin < 3
    DataStr = [];
end

% Set channel locations from template
[chanlocs, ndchanlocs] = template_to_chanlocs(which('GSN-HydroCel-257.sfp'));
[~, system] = ishdeeg({chanlocs.labels});
incl = hdeeg_scalpchannels(system);
idx_keep = ismember({chanlocs.labels}, incl);
chanlocs = chanlocs(idx_keep);

% Extract the Key-Value pairs from the filename
kv = filename2struct(reverse_fileparts(outfname));
% Init empty first-level output structure
EEG = emptyfstlvl();
% Set meta info
EEG.filename = [struct2filename(kv), '.mat'];
EEG.filepath = sprintf('%s/derivatives/EEG-output-fstlvl/sub-%s/ses-%s', strrep(pwd, '\', '/'), kv.sub, kv.ses);
EEG.subject = kv.sub;
EEG.session = kv.ses;
EEG.task = kv.task;
EEG.run = '1';
EEG.condition = kv.ses;
% Set data info
EEG.nbchan = length(chanlocs);
EEG.srate = 0;
% Set Data
if ~isempty(DataStr)
    flds = fieldnames(DataStr);
    for i = 1:length(flds)
        EEG.(flds{i}) = DataStr.(flds{i});
    end
end
% Set features (DV for group analysis)
EEG.features = Features;
% Check channels
for i = 1:length(EEG.features)
    if length(EEG.features(i).data) ~= 178
        error('Channels not right')
    end
end
% Set channel locations from template
EEG.chanlocs = chanlocs;
EEG.chaninfo.ndchanlocs = ndchanlocs;
% Set JSON info
EEG.etc.JSON.Description = 'Abc';
EEG.etc.JSON.Sources = 'unknown';
EEG.etc.JSON.TaskName = kv.task;
EEG.etc.JSON.EEGReference = 'average';
EEG.etc.JSON.EEGChannelCount = EEG.nbchan;
EEG.etc.JSON.ECGChannelCount = 0;
EEG.etc.JSON.EMGChannelCount = 0;
EEG.etc.JSON.EOGChannelCount = 0;
EEG.etc.JSON.MiscChannelCount = 0;
EEG.etc.JSON.TrialCount = 1;
% Save
if nargout == 0
    SaveDataset(EEG, 'matrix');
end

end