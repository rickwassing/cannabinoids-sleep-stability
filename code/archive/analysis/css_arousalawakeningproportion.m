clear
clc
ARO = readtable('phenotype/2024-07-26T105951_arousals.csv');

ARO.stage_str = cell(size(ARO, 1), 1);
ARO.stage_str(ARO.stage == 1) = {'wake'};
ARO.stage_str(ARO.stage == 0) = {'rem'};
ARO.stage_str(ARO.stage == -1) = {'n1'};
ARO.stage_str(ARO.stage == -2) = {'n2'};
ARO.stage_str(ARO.stage == -3) = {'n3'};

idx = ARO.stage < 0;
mdl_nrem = fitglme(ARO(idx, :), 'is_awakening ~ 1 + aro_type + cond + (1 | participant_id)', 'Distribution', 'binomial')

idx = ARO.stage == 0;
mdl_rem = fitglme(ARO(idx, :), 'is_awakening ~ 1 + aro_type + cond + (1 | participant_id)', 'Distribution', 'binomial')