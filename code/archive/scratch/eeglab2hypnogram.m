% ARCHIVED 2026-09-18: byte-for-byte duplicate of utils/eeg-io/css_eeglab2hypnogram.m
% (SRP_REFACTOR_PLAN.md Phase 2). Its only caller,
% subject-level/spectral-features/css_extractspindles.m, was pointed at
% css_eeglab2hypnogram instead (the version every other caller in code/ already used).
% Kept for provenance, not called from main.m.
% -------------------------------------------------------------------------
function HYP = eeglab2hypnogram(EEG)
% Extract sleep stages
HYP = EEG.event(ismember(lower({EEG.event.type}), {'n1', 'n2', 'n3', 'w', 'r', 'wake', 'rem', 's1', 's2', 's3', 'nrem1', 'nrem2', 'nrem3', '1', '2', '3'}));
stages = {HYP.type};
% Calculte sleep cycles
HYP = sleep_cycles(stages, {'wake', 'n1', 'n2', 'n3', 'rem'});
% Transform data and create output table
HYP.times = HYP.times';
HYP.stage_num = HYP.sleepstage'; 
HYP.stage_str = stages';
HYP.cycle = HYP.sleepcycle';
HYP.episode = HYP.sleepepisode';
HYP.latency = HYP.times * 24*60*60*EEG.srate;
HYP.sleepwin = false(length(HYP.times), 1);
HYP.sleepwin(find(HYP.sleepstage ~= 1, 1, 'first') : find(HYP.sleepstage ~= 1, 1, 'last')) = true;
HYP = rmfield(HYP, {'sleepcycle', 'sleepepisode', 'sleepstage'});
HYP = struct2table(HYP);
end