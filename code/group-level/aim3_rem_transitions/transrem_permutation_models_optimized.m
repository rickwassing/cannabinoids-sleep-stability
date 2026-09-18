function [perms] = transrem_permutation_models_optimized(T, tme)

nperms = 100;

% ------------------------------------------------------------------------
% Build base design table (adding y and Pr columns) if not already present
% -------------------------------------------------------------------------
desmat = table();
desmat.cond = T.cond;
desmat.sub = T.sub;

idx_t = find(tme >= -90 & tme < 0);
num_timepoints = length(idx_t);

% Precompute y and Pr columns into the table only once:
for tme = 1:num_timepoints
    desmat.(sprintf('y%i', tme))  = cellfun(@(d) double(d(tme)), T.s);
end

desmat_orig = desmat;

% ------------------------------------------------------------------------
% Precompute model spec strings
% -------------------------------------------------------------------------
mspecs.y.cond  = arrayfun(@(t) sprintf('y%i ~ 1 + cond + (1|sub)', t), 1:num_timepoints, 'UniformOutput', false);

% ------------------------------------------------------------------------
% Fit ONE reference model to extract design matrices
% ------------------------------------------------------------------------
mdl_ref = fitlme(desmat, mspecs.y.cond{1});

X = mdl_ref.designMatrix('Fixed');
Z = X(:, 1); % random intercept
group = nominal(desmat.sub);

% ------------------------------------------------------------------------
% Preallocate perms struct array (single template then repmat)
% -------------------------------------------------------------------------
% create a template for one permutation result
emptyPerm.y.tstat  = nan(1, num_timepoints);
emptyPerm.y.pval  = nan(1, num_timepoints);

perms = repmat(emptyPerm, nperms, 1);

% -------------------------------------------------------------------------
% Main permutation loop
% -------------------------------------------------------------------------
rt = now; %#ok<TNOW1>
for p = 1:nperms

    % Start each permutation with original data
    desmat = desmat_orig;

    % For p>1 perform label permutation(s)
    if p > 1
        desmat = permute_event_labels(desmat, 'cond');
    end

    for tme = 1:num_timepoints
        y = desmat.(sprintf('y%i', tme));
        try
            %mdl = fitlme(desmat, mspecs.y.cond{tme});
            mdl = fitlmematrix(X, y, Z, group, 'FitMethod', 'ML', 'CovariancePattern', 'FullCholesky');
            perms(p).y.tstat(tme) = mdl.Coefficients.tStat(2);
            perms(p).y.pval(tme) = mdl.Coefficients.pValue(2);
        catch ME
        end
        
        rt = remainingTime(rt, nperms*num_timepoints, true);

    end % time loop

end % perm loop

end
