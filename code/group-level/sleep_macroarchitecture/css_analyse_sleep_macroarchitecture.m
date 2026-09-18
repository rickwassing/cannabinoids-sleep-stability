function css_analyse_sleep_macroarchitecture()
% Load data
PSG = readtable('phenotype/2024-07-26T171508_psg-variables.csv');
ARO = readtable('phenotype/2024-07-26T171508_arousals.csv');
ARO.duration(ARO.duration > 30) = 30;  % Cap extreme durations

% Stage coding
ARO.stage_str = cell(size(ARO,1),1);
ARO.stage_str(ARO.stage ==  1) = {'r'};
ARO.stage_str(ARO.stage ==  0) = {'w'};
ARO.stage_str(ARO.stage == -1) = {'n1'};
ARO.stage_str(ARO.stage == -2) = {'n2'};
ARO.stage_str(ARO.stage == -3) = {'n3'};

% ------------------------------------------------------------------------
%  Compare PSG variables between conditions (Table 1)
% -------------------------------------------------------------------------
compare_psg_variables(PSG);

%% ------------------------------------------------------------------------
%  Mixed-effects modelling: Onset, Duration, Awakening probability
% -------------------------------------------------------------------------
[ARO, S] = model_awakening_probability(ARO);

%% ------------------------------------------------------------------------
%  Figure: Arousal descriptives
% -------------------------------------------------------------------------
plot_arousal_descriptives(PSG, ARO, S);

end
