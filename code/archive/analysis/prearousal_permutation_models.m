function [perms] = prearousal_permutation_models(T)

clear mspecs

perms = struct();

desmat = table();
desmat.cond = T.ses;
desmat.is_awakening = T.is_awakening;
desmat.sub = T.sub;
desmat.chan = T.channel;
desmat.id = T.id;
desmat.y_amp = T.amp_d0;

for t = 1:600
    desmat.(sprintf('y%i', t)) = cellfun(@(d) double(d(t)), T.smtdata);
    desmat.(sprintf('Pr%i', t)) = cellfun(@(d) double(d(t)), T.pr_aro);
end
desmat = averageSigmaPerEvent(desmat);
desmat_orig = desmat;

mspecs.y.cond = arrayfun(@(t) sprintf('y%i ~ 1 + cond + (1|sub)', t), 1:600, 'UniformOutput', false);
mspecs.y.awake = arrayfun(@(t) sprintf('y%i ~ 1 + is_awakening + (1|sub)', t), 1:600, 'UniformOutput', false);
mspecs.pr.cond = arrayfun(@(t) sprintf('Pr%i ~ 1 + cond + (1|sub)', t), 1:600, 'UniformOutput', false);
mspecs.pr.awake = arrayfun(@(t) sprintf('Pr%i ~ 1 + is_awakening + (1|sub)', t), 1:600, 'UniformOutput', false);

rt = now; %#ok<TNOW1>
nperms = 100;

for p = 1:nperms

    desmat = desmat_orig;

    if p > 1
        desmat = permuteEventLabels(desmat, 'cond');
        desmat = permuteEventLabels(desmat, 'is_awakening');
    end

    perms(p).y.aw.tstat = nan(1, 600);
    perms(p).y.aw.pval = nan(1, 600);
    perms(p).y.cs.tstat = nan(1, 600);
    perms(p).y.cs.pval = nan(1, 600);
    perms(p).y.pbo.tstat = nan(1, 600);
    perms(p).y.pbo.pval = nan(1, 600);
    perms(p).y.etc.tstat = nan(1, 600);
    perms(p).y.etc.pval = nan(1, 600);
    
    perms(p).pr.aw.tstat = nan(1, 600);
    perms(p).pr.aw.pval = nan(1, 600);
    perms(p).pr.cs.tstat = nan(1, 600);
    perms(p).pr.cs.pval = nan(1, 600);
    perms(p).pr.pbo.tstat = nan(1, 600);
    perms(p).pr.pbo.pval = nan(1, 600);
    perms(p).pr.etc.tstat = nan(1, 600);
    perms(p).pr.etc.pval = nan(1, 600);

    desmat_aw_true = desmat(desmat.is_awakening == "true", :);
    desmat_aw_false = desmat(desmat.is_awakening == "false", :);

    for t = 1:600

        try
            % Test differences in sigma between PBO and ETC in AW cond
            mdl = fitlme(desmat_aw_true, mspecs.y.cond{t});
            perms(p).y.aw.tstat(t) = mdl.Coefficients.tStat(2);
            perms(p).y.aw.pval(t) = mdl.Coefficients.pValue(2);
            % Test differences in sigma between PBO and ETC in CS cond
            mdl = fitlme(desmat_aw_false, mspecs.y.cond{t});
            perms(p).y.cs.tstat(t) = mdl.Coefficients.tStat(2);
            perms(p).y.cs.pval(t) = mdl.Coefficients.pValue(2);
        catch ME
        end

        % try
        %     % Test differences in sigma between AW and CS in PBO
        %     mdl = fitglme(desmat(desmat.cond == "placebo", :), mspecs.y.awake{t}, 'Distribution', 'Normal');
        %     perms(p).y.pbo.tstat(t) = mdl.Coefficients.tStat(2);
        %     perms(p).y.pbo.pval(t) = mdl.Coefficients.pValue(2);
        %     % Test differences in sigma between AW and CS in ETC
        %     mdl = fitglme(desmat(desmat.cond == "etc120", :), mspecs.y.awake{t}, 'Distribution', 'Normal');
        %     perms(p).y.etc.tstat(t) = mdl.Coefficients.tStat(2);
        %     perms(p).y.etc.pval(t) = mdl.Coefficients.pValue(2);
        % catch ME
        % end

        try
            % Test differences in Pr(aro) between PBO and ETC in AW cond
            mdl = fitglme(desmat_aw_true, mspecs.pr.cond{t}, 'Distribution', 'Binomial');
            perms(p).pr.aw.tstat(t) = mdl.Coefficients.tStat(2);
            perms(p).pr.aw.pval(t) = mdl.Coefficients.pValue(2);
            % Test differences in Pr(aro) between PBO and ETC in CS cond
            mdl = fitglme(desmat_aw_false, mspecs.pr.cond{t}, 'Distribution', 'Binomial');
            perms(p).pr.cs.tstat(t) = mdl.Coefficients.tStat(2);
            perms(p).pr.cs.pval(t) = mdl.Coefficients.pValue(2);
        catch ME
        end

        % try
        %     % Test differences in Pr(aro) between AW and CS in PBO
        %     mdl = fitglme(desmat(desmat.cond == "placebo", :), mspecs.pr.awake{t}, 'Distribution', 'Binomial');
        %     perms(p).pr.pbo.tstat(t) = mdl.Coefficients.tStat(2);
        %     perms(p).pr.pbo.pval(t) = mdl.Coefficients.pValue(2);
        %     % Test differences in Pr(aro) between AW and CS in ETC
        %     mdl = fitglme(desmat(desmat.cond == "etc120", :), mspecs.pr.awake{t}, 'Distribution', 'Binomial');
        %     perms(p).pr.etc.tstat(t) = mdl.Coefficients.tStat(2);
        %     perms(p).pr.etc.pval(t) = mdl.Coefficients.pValue(2);
        % catch ME
        % end

        rt = remainingTime(rt, nperms*600, true);

    end


end

end