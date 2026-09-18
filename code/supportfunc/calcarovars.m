function t = calcarovars(psgvars, aro)

t = table();

t.n_aro = sum(strcmpi(aro.aro_type, '1arousal'));
t.n_alp = sum(strcmpi(aro.aro_type, '2alpha'));
t.n_mcr = sum(strcmpi(aro.aro_type, '3micro'));

t.n_aro_nrem = sum(strcmpi(aro.aro_type(aro.stage < 0), '1arousal'));
t.n_aro_n1 = sum(strcmpi(aro.aro_type(aro.stage == -1), '1arousal'));
t.n_aro_n2 = sum(strcmpi(aro.aro_type(aro.stage == -2), '1arousal'));
t.n_aro_n3 = sum(strcmpi(aro.aro_type(aro.stage == -3), '1arousal'));
t.n_aro_rem = sum(strcmpi(aro.aro_type(aro.stage == 0), '1arousal'));

t.n_alp_nrem = sum(strcmpi(aro.aro_type(aro.stage < 0), '2alpha'));
t.n_alp_n1 = sum(strcmpi(aro.aro_type(aro.stage == -1), '2alpha'));
t.n_alp_n2 = sum(strcmpi(aro.aro_type(aro.stage == -2), '2alpha'));
t.n_alp_n3 = sum(strcmpi(aro.aro_type(aro.stage == -3), '2alpha'));
t.n_alp_rem = sum(strcmpi(aro.aro_type(aro.stage == 0), '2alpha'));

t.n_mcr_nrem = sum(strcmpi(aro.aro_type(aro.stage < 0), '3micro'));
t.n_mcr_n1 = sum(strcmpi(aro.aro_type(aro.stage == -1), '3micro'));
t.n_mcr_n2 = sum(strcmpi(aro.aro_type(aro.stage == -2), '3micro'));
t.n_mcr_n3 = sum(strcmpi(aro.aro_type(aro.stage == -3), '3micro'));
t.n_mcr_rem = sum(strcmpi(aro.aro_type(aro.stage == 0), '3micro'));

t.i_aro = 60 * t.n_aro / psgvars.tst;
t.i_alp = 60 * t.n_alp / psgvars.tst;
t.i_mcr = 60 * t.n_mcr / psgvars.tst;

t.i_aro_nrem = 60 * t.n_aro_nrem / (psgvars.n1_min + psgvars.n2_min + psgvars.n3_min);
t.i_aro_n1 = 60 * t.n_aro_n1 / psgvars.n1_min;
t.i_aro_n2 = 60 * t.n_aro_n2 / psgvars.n2_min;
t.i_aro_n3 = 60 * t.n_aro_n3 / psgvars.n3_min;
t.i_aro_rem = 60 * t.n_aro_rem / psgvars.rem_min;

t.i_alp_nrem = 60 * t.n_alp_nrem / (psgvars.n1_min + psgvars.n2_min + psgvars.n3_min);
t.i_alp_n1 = 60 * t.n_alp_n1 / psgvars.n1_min;
t.i_alp_n2 = 60 * t.n_alp_n2 / psgvars.n2_min;
t.i_alp_n3 = 60 * t.n_alp_n3 / psgvars.n3_min;
t.i_alp_rem = 60 * t.n_alp_rem / psgvars.rem_min;

t.i_mcr_nrem = 60 * t.n_mcr_nrem / (psgvars.n1_min + psgvars.n2_min + psgvars.n3_min);
t.i_mcr_n1 = 60 * t.n_mcr_n1 / psgvars.n1_min;
t.i_mcr_n2 = 60 * t.n_mcr_n2 / psgvars.n2_min;
t.i_mcr_n3 = 60 * t.n_mcr_n3 / psgvars.n3_min;
t.i_mcr_rem = 60 * t.n_mcr_rem / psgvars.rem_min;

end