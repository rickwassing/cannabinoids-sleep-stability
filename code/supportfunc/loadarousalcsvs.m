function ARO = loadarousalcsvs(COND, varargin)
selectcondition = false;
for i = 1:2:length(varargin)
    switch lower(varargin{i})
        case 'selectcondition'
            selectcondition = varargin{i+1};
    end
end

Arousals = dir('derivatives/EEG-preproc/sub-*/ses-*/sub*_arousal.csv');
ARO = [];
for i = 1:length(Arousals)
    kv = filename2struct(Arousals(i).name);
    cond = COND.condition{strcmpi(COND.folder_name, kv.sub)};
    if strcmpi(cond, 'A') && selectcondition % Condition B is placebo
        continue
    end
    a = readtable(fullfile(Arousals(i).folder, Arousals(i).name));
    jfile = dir([Arousals(i).folder, '/sub*_eeg.json']);
    j = json2struct(fullfile(jfile(1).folder, jfile(1).name));
    a.latency = a.latency ./ j.SamplingFrequency;
    a.duration = a.duration ./ j.SamplingFrequency;
    a.interval = a.interval ./ j.SamplingFrequency;
    a.filename = repmat({Arousals(i).name}, size(a, 1), 1);
    a = movevars(a, 'filename', 'Before', 'latency');
    if isempty(ARO)
        ARO = a;
    else
        ARO = [ARO; a];
    end
end
% rmidx = ARO.interval < 1;
% ARO(rmidx, :) = [];


% Remove entries we know are missing from the ARO epoched files
rmIdx = ...
    contains(ARO.filename, 'sub-r0122lg') & ARO.id == 1605;
ARO(rmIdx, :) = [];

ARO.newtype = ARO.type;
ARO.newtype(strcmpi(ARO.newtype, 'arousalemg') & ARO.duration < 3) = {'microarousalemg'};
ARO.newtype(strcmpi(ARO.newtype, 'arousal') & ARO.duration < 3) = {'microarousal'};
ARO.newtype(strcmpi(ARO.newtype, 'alphaemg') & ARO.duration < 3) = {'microalphaemg'};
ARO.newtype(strcmpi(ARO.newtype, 'alpha') & ARO.duration < 3) = {'microalpha'};
ARO.inttype = ARO.newtype;
ARO.inttype(contains(ARO.inttype, 'micro')) = {'3micro'};
ARO.inttype(contains(ARO.inttype, 'alpha')) = {'2alpha'};
ARO.inttype(contains(ARO.inttype, 'arousal')) = {'1arousal'};
ARO.emgtype = repmat({'2no-emg'}, size(ARO, 1), 1);
ARO.emgtype(contains(ARO.newtype, 'emg')) = {'1emg'};
end