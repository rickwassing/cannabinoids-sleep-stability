function css_analyse_1a()
% -------------------------------------------------------------------------
% STEP 1:
% Use 300-second bouts of continuous N2 sleep to describe the sigma/spindle
% ISF between conditions:
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% SUB: a
% - Describe selected bouts (average number, duration, number of interrupting
%   epochs of each sleep stage.
clc
Files = dir('derivatives/EEG-segmented/sub-*/ses-*/sub-*_desc-sigmanrembout*.set');
TBL = table();
for f = 1:length(Files)
    % Load header file
    EEG = LoadDataset(fullfile(Files(f).folder, Files(f).name), 'header');
    kv = filename2struct(EEG.setname);
    % Init table
    tbl = table();
    tbl.sub = {kv.sub};
    tbl.ses = {kv.ses};
    tbl.fname = {EEG.setname};
    % Find bouts
    idx_b = [ ...
        find(strcmpi({EEG.event.type}, 'n2'), 1, 'first'), ...
        find(strcmpi({EEG.event.type}, 'boundary')), ...
        find(strcmpi({EEG.event.type}, 'n2'), 1, 'last')];
    tbl.nbouts = length(idx_b)-1;
    % Extract duration, number of interrupting epochs
    dur = [];
    int = struct();
    int.wake = 0;
    int.n1 = 0;
    int.n3 = 0;
    int.rem = 0;
    for b = 1:length(idx_b)-1
        latency = [EEG.event(idx_b(b):idx_b(b+1)).origlatency];
        stages = {EEG.event(idx_b(b):idx_b(b+1)).type};
        idx = ismember(stages, {'n1', 'n2', 'n3', 'rem', 'wake'});
        latency = latency(idx);
        stages = stages(idx);
        idx = find(strcmpi(stages, 'n2'), 1, 'first'):find(strcmpi(stages, 'n2'), 1, 'last');
        latency = latency(idx);
        stages = stages(idx);
        dur(b) = (latency(end) - latency(1)) / EEG.srate + 30;
        int.wake = int.wake + sum(diff([0, find(strcmpi(stages, 'wake'))]) > 1);
        int.n1 = int.n1 + sum(diff([0, find(strcmpi(stages, 'n1'))]) > 1);
        int.n3 = int.n3 + sum(diff([0, find(strcmpi(stages, 'n3'))]) > 1);
        int.rem = int.rem + sum(diff([0, find(strcmpi(stages, 'rem'))]) > 1);
    end
    % Average duration
    tbl.av_dur = mean(dur);
    % Average number of interrupting epochs
    tbl.intr_wake = int.wake / (length(idx_b)-1);
    tbl.intr_n1 = int.n1 / (length(idx_b)-1);
    tbl.intr_n3 = int.n3 / (length(idx_b)-1);
    tbl.intr_rem = int.rem / (length(idx_b)-1);
    % Assign to table
    if f == 1
        TBL = tbl;
    else
        TBL = [TBL; tbl];
    end
end
mdl.nbouts = fitglme(TBL, 'nbouts ~ 1 + ses + (1 | sub)', 'Distribution', 'Poisson');
mdl.duration = fitglme(TBL, 'av_dur ~ 1 + ses + (1 | sub)');
fprintf('##############################################################\n')
fprintf('# REPORT\n')
fprintf('Continuous bouts of NREM stage-2 sleep were selected\n')
fprintf('that contained consecutive stage-2 epochs and were interrupted\n')
fprintf('by at most one epoch of any other sleep stage. On average %.1f (%.1f) bouts were\n', mean(TBL.nbouts), std(TBL.nbouts))
fprintf('selected from each recording with a duration of %.0f (%.0f) seconds. There were no significant differences\n', mean(TBL.av_dur), std(TBL.av_dur))
fprintf('between conditions in the number of selected bouts (t(%i) = %.1f, p = %.2f),\n', mdl.nbouts.Coefficients.DF(2), mdl.nbouts.Coefficients.tStat(2), mdl.nbouts.Coefficients.pValue(2))
fprintf('or their duration (t(%i) = %.1f, p = %.2f).\n', mdl.duration.Coefficients.DF(2), mdl.duration.Coefficients.tStat(2), mdl.duration.Coefficients.pValue(2))
fprintf('On average, the bouts were interrupted by %.1f (%.1f) wake epochs,\n', mean(TBL.intr_wake), std(TBL.intr_wake))
fprintf('%.1f (%.1f) N1 epochs,\n', mean(TBL.intr_n1), std(TBL.intr_n1))
fprintf('%.1f (%.1f) N3 epochs, and\n', mean(TBL.intr_n3), std(TBL.intr_n3))
fprintf('%.1f (%.1f) REM epochs.\n', mean(TBL.intr_rem), std(TBL.intr_rem))
fprintf('##############################################################\n')

end