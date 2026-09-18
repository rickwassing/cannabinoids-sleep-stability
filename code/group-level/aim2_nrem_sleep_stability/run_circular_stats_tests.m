function Circ = run_circular_stats_tests(T_this, pcfg, roi, delay, Circ)
% -------------------------------------------------------------------------
% Run the circular-statistics test battery (Rayleigh, Hermans-Rasson, Dip,
% peak-density, Kuiper, Watson's U2) comparing phase angles between
% awakening/continued-sleep arousals and placebo/THC-CBD conditions for a
% single region of interest ('roi'), and store the parsed p-values in the
% running 'Circ' table.
% -------------------------------------------------------------------------
r = 0;
for fld = {'aw', 'cs'}

    fprintf('---------------------\n');
    fprintf('%s\n', fld{:});

    for cond = {'pbo', 'etc'}
        r = r+1;
        AData = within_chan_circ_mean(T_this(pcfg.idx.(cond{:}).(fld{:}), :), delay);

        [pval, m] = circ_rtest(AData);
        fprintf('Rayleigh test for non-uniformity of ''%s'' arousals in ''%s'' condition (m = %.2f, p = %.3f).\n', fld{:}, cond{:}, m, pval)
        Circ.(['Ray_', roi]){r} = parsepvalue(pval);

        [pval, T] = circ_hrtest(AData);
        fprintf('Hermans-Rasson test for non-uniformity of ''%s'' arousals in ''%s'' condition (m = %.2f, p = %.3f).\n', fld{:}, cond{:}, T, pval)
        Circ.(['HR_', roi]){r} = parsepvalue(pval);

        [pval, dip, xl, xu] = circ_diptest(AData);
        fprintf('Dip-test of ''%s'' arousals in ''%s'' condition (dip = %.2f, p = %.3f, limits %.2f - %.2f).\n', fld{:}, cond{:}, dip, pval, circ_rad2deg360(xl), circ_rad2deg360(xu))
        Circ.(['Dip_', roi]){r} = parsepvalue(pval);

        WData = -pi:pi/180:pi;
        KData = circ_ksdensity(AData, WData, [-pi, pi]);
        [pks, plocs, pwidth] = findpeaks(KData, 'SortStr','descend'); %#ok<ASGLU>
        if pval > 0.05
            fprintf('Unimodal peak at %.2f (%.2f - %.2f) degrees\n', circ_rad2deg360(WData(plocs(1))), circ_rad2deg360(WData(plocs(1))-pwidth(1)*mean(diff(WData))), circ_rad2deg360(WData(plocs(1))+pwidth(1)*mean(diff(WData))))
        else
            fprintf('Multimodal peaks at %.2f (%.2f - %.2f) and %.2f (%.2f - %.2f) degrees\n', ...
                circ_rad2deg360(WData(plocs(1))), circ_rad2deg360(WData(plocs(1))-pwidth(1)*mean(diff(WData))), circ_rad2deg360(WData(plocs(1))+pwidth(1)*mean(diff(WData))), ...
                circ_rad2deg360(WData(plocs(2))), circ_rad2deg360(WData(plocs(2))-pwidth(2)*mean(diff(WData))), circ_rad2deg360(WData(plocs(2))+pwidth(2)*mean(diff(WData))))
        end

        fprintf('\n')
    end
end

for fld = {'aw', 'cs'}
    [pval, k] = circ_kuipertest(...
        within_chan_circ_mean(T_this(pcfg.idx.pbo.(fld{:}), :), delay), ...
        within_chan_circ_mean(T_this(pcfg.idx.etc.(fld{:}), :), delay), ...
        pcfg.nbins, false);
    fprintf('Kuiper-test indicated phase angles are different between PBO and ETC for ''%s'' arousals (k = %.2f, p = %.3f).\n', fld{:}, k, pval)

    A1 = within_chan_circ_mean(T_this(pcfg.idx.pbo.(fld{:}), :), delay);
    A2 = within_chan_circ_mean(T_this(pcfg.idx.etc.(fld{:}), :), delay);
    [pval, U2_obs, U2_H0] = watsons_U2_perm_test(A1,A2, 200); %#ok<ASGLU>
    fprintf('Nonparametric permutation test based on Watson''s U2 indicated phase angles are/are not different between PBO and ETC for ''%s'' arousals (U2 = %.2f, p = %.3f).\n', fld{:}, U2_obs, pval)

end

end
