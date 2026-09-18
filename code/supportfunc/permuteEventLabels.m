function T = permuteEventLabels(T, labelVarName)
%PERMUTEEVENTLABELS Permute 'cond' or 'is_awakening' respecting subject structure
%
% INPUTS:
%   T            - event-level table (one row per arousal)
%   labelVarName - 'cond' or 'is_awakening'
%
% OUTPUT:
%   T            - table with new column [labelVarName '_perm']

% Ensure categorical
T.sub     = categorical(T.sub);
T.cond    = categorical(T.cond);
try
T.is_awakening = categorical(T.is_awakening);
end

switch labelVarName
    case 'is_awakening'
        % Shuffle is_awakening within each subject
        label_perm = T.is_awakening; % initialize
        subjects = categories(T.sub);
        for i = 1:numel(subjects)
            idx = find(T.sub == subjects{i});
            label_perm(idx) = label_perm(idx(randperm(length(idx))));
        end

    case 'cond'
        % Paired swap per subject (crossover design)
        label_perm = T.cond; % initialize
        subjects = categories(T.sub);
        for i = 1:numel(subjects)
            idx = T.sub == subjects{i};
            if rand > 0.5
                % swap drug <-> placebo
                tmp = label_perm(idx);
                tmp(tmp=="etc120") = "tmp";
                tmp(tmp=="placebo") = "etc120";
                tmp(tmp=="tmp") = "placebo";
                label_perm(idx) = tmp;
            end
        end

    otherwise
        error('labelVarName must be ''cond'' or ''is_awakening''');
end

% Add permuted column to table
T.(labelVarName) = label_perm;

end