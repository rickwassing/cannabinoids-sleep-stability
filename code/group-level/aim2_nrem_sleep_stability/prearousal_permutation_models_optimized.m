function [perms] = prearousal_permutation_models_optimized(T)

nperms = 25;

% ------------------------------------------------------------------------
% Build base design table (adding y and Pr columns) if not already present
% -------------------------------------------------------------------------
desmat = table();
desmat.cond = T.ses;
desmat.is_awakening = T.is_awakening;
desmat.sub = T.sub;
desmat.chan = T.channel;
desmat.id = T.id;
desmat.y_amp = T.amp_d0;

num_timepoints = 600;

% Precompute y and Pr columns into the table only once:
for t = 1:num_timepoints
    desmat.(sprintf('y%i', t))  = cellfun(@(d) double(d(t)), T.smtdata);
    desmat.(sprintf('Pr%i', t)) = cellfun(@(d) double(d(t)), T.pr_aro);
end

desmat = average_sigma_per_event(desmat);
desmat_orig = desmat;

% ------------------------------------------------------------------------
% Precompute model spec strings
% -------------------------------------------------------------------------
mspecs.y.cond  = arrayfun(@(t) sprintf('y%i ~ 1 + cond + (1|sub)', t), 1:num_timepoints, 'UniformOutput', false);
mspecs.y.awake = arrayfun(@(t) sprintf('y%i ~ 1 + is_awakening + (1|sub)', t), 1:num_timepoints, 'UniformOutput', false);
mspecs.pr.cond  = arrayfun(@(t) sprintf('Pr%i ~ 1 + cond + (1|sub)', t), 1:num_timepoints, 'UniformOutput', false);
mspecs.pr.awake = arrayfun(@(t) sprintf('Pr%i ~ 1 + is_awakening + (1|sub)', t), 1:num_timepoints, 'UniformOutput', false);

% ------------------------------------------------------------------------
% Preallocate perms struct array (single template then repmat)
% -------------------------------------------------------------------------
% create a template for one permutation result
emptyPerm.y.aw.tstat  = nan(1, num_timepoints);
emptyPerm.y.aw.pval  = nan(1, num_timepoints);
emptyPerm.y.cs.tstat  = nan(1, num_timepoints);
emptyPerm.y.cs.pval  = nan(1, num_timepoints);
emptyPerm.y.pbo.tstat = nan(1, num_timepoints);
emptyPerm.y.pbo.pval = nan(1, num_timepoints);
emptyPerm.y.etc.tstat = nan(1, num_timepoints);
emptyPerm.y.etc.pval = nan(1, num_timepoints);

emptyPerm.pr.aw.tstat  = nan(1, num_timepoints);
emptyPerm.pr.aw.pval  = nan(1, num_timepoints);
emptyPerm.pr.cs.tstat  = nan(1, num_timepoints);
emptyPerm.pr.cs.pval  = nan(1, num_timepoints);
emptyPerm.pr.pbo.tstat = nan(1, num_timepoints);
emptyPerm.pr.pbo.pval = nan(1, num_timepoints);
emptyPerm.pr.etc.tstat = nan(1, num_timepoints);
emptyPerm.pr.etc.pval = nan(1, num_timepoints);

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
        desmat = permute_event_labels(desmat, 'is_awakening');
    end

    desmat_aw_true = desmat(desmat.is_awakening == "true", :);
    desmat_aw_false = desmat(desmat.is_awakening == "false", :);

    % If there are very small groups, subsequent fits may fail; still we proceed
    for t = 1:num_timepoints
        % ---------- Y models (LME via fitlme) ----------
        try
            % AW: PBO vs ETC in awakening condition
            if p > 1
                if any(perms(1).y.aw.pval < 0.05)
                    mdl = fitlme(desmat_aw_true, mspecs.y.cond{t});
                    perms(p).y.aw.tstat(t) = mdl.Coefficients.tStat(2);
                    perms(p).y.aw.pval(t) = mdl.Coefficients.pValue(2);
                end
            else
                mdl = fitlme(desmat_aw_true, mspecs.y.cond{t});
                perms(p).y.aw.tstat(t) = mdl.Coefficients.tStat(2);
                perms(p).y.aw.pval(t) = mdl.Coefficients.pValue(2);
            end
            % CS: PBO vs ETC in not-awakening condition
            if p > 1
                if any(perms(1).y.cs.pval < 0.05)
                    mdl = fitlme(desmat_aw_false, mspecs.y.cond{t});
                    perms(p).y.cs.tstat(t) = mdl.Coefficients.tStat(2);
                    perms(p).y.cs.pval(t) = mdl.Coefficients.pValue(2);
                end
            else
                mdl = fitlme(desmat_aw_false, mspecs.y.cond{t});
                perms(p).y.cs.tstat(t) = mdl.Coefficients.tStat(2);
                perms(p).y.cs.pval(t) = mdl.Coefficients.pValue(2);
            end
        catch ME
        end

        % ---------- Pr models (Binomial GLME via fitglme) ----------
        try
            % AW condition
            if p > 1
                if any(perms(1).pr.aw.pval < 0.05)
                    mdl = fitglme(desmat_aw_true, mspecs.pr.cond{t}, 'Distribution', 'Binomial');
                    perms(p).pr.aw.tstat(t) = mdl.Coefficients.tStat(2);
                    perms(p).pr.aw.pval(t)  = mdl.Coefficients.pValue(2);
                end
            else
                mdl = fitglme(desmat_aw_true, mspecs.pr.cond{t}, 'Distribution', 'Binomial');
                perms(p).pr.aw.tstat(t) = mdl.Coefficients.tStat(2);
                perms(p).pr.aw.pval(t)  = mdl.Coefficients.pValue(2);
            end
            % CS condition
            if p > 1
                if any(perms(1).pr.cs.pval < 0.05)
                    mdl = fitglme(desmat_aw_false, mspecs.pr.cond{t}, 'Distribution', 'Binomial');
                    perms(p).pr.cs.tstat(t) = mdl.Coefficients.tStat(2);
                    perms(p).pr.cs.pval(t)  = mdl.Coefficients.pValue(2);
                end
            else
                mdl = fitglme(desmat_aw_false, mspecs.pr.cond{t}, 'Distribution', 'Binomial');
                perms(p).pr.cs.tstat(t) = mdl.Coefficients.tStat(2);
                perms(p).pr.cs.pval(t)  = mdl.Coefficients.pValue(2);
            end
        catch ME
            % keep NaNs if model fails
        end

        rt = remainingTime(rt, nperms*num_timepoints, true);

    end % time loop

end % perm loop

end
