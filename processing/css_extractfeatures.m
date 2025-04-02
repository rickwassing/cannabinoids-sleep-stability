function features = css_extractfeatures(FIT, ISF)

features = struct();
i = 0;

i = i+1;
features(i).label = 'isffreq';
features(i).type = 'hz';
features(i).data = ascolumn([FIT.mu]);

i = i+1;
features(i).label = 'isfpeak';
features(i).type = 'au';
features(i).data = ascolumn([FIT.peak]);

i = i+1;
features(i).label = 'isfbandwidth';
features(i).type = 'hz';
features(i).data = ascolumn([FIT.bandwidth]);

i = i+1;
features(i).label = 'isftotpow';
features(i).type = 'au';
features(i).data = ascolumn(sum(cat(2, FIT.pow), 1));

i = i+1;
idx_f = ISF.freqs > 0.0099 & ISF.freqs < 0.0268;
features(i).label = sprintf('%.1f-%.1fs', min(1./ISF.freqs(idx_f)), max(1./ISF.freqs(idx_f)));
features(i).type = 'au';
features(i).data = mean(ISF.data(:, idx_f), 2);

i = i+1;
idx_f = ISF.freqs > 0.0266 & ISF.freqs < 0.0401;
features(i).label = sprintf('%.1f-%.1fs', min(1./ISF.freqs(idx_f)), max(1./ISF.freqs(idx_f)));
features(i).type = 'au';
features(i).data = mean(ISF.data(:, idx_f), 2);

i = i+1;
idx_f = ISF.freqs > 0.0399 & ISF.freqs < 0.0534;
features(i).label = sprintf('%.1f-%.1fs', min(1./ISF.freqs(idx_f)), max(1./ISF.freqs(idx_f)));
features(i).type = 'au';
features(i).data = mean(ISF.data(:, idx_f), 2);

i = i+1;
idx_f = ISF.freqs > 0.0532 & ISF.freqs < 0.1001;
features(i).label = sprintf('%.1f-%.1fs', min(1./ISF.freqs(idx_f)), max(1./ISF.freqs(idx_f)));
features(i).type = 'au';
features(i).data = mean(ISF.data(:, idx_f), 2);

end