function ISF = interp_isfspectrum(ISF, badchans)

for i = 1:length(badchans)
    p = [ISF.chanlocs(badchans(i)).X; ISF.chanlocs(badchans(i)).Y; ISF.chanlocs(badchans(i)).Z];
    q = [[ISF.chanlocs.X]; [ISF.chanlocs.Y]; [ISF.chanlocs.Z]];
    d = ((p(1)-q(1, :)).^2 + (p(2)-q(2, :)).^2 + (p(3)-q(3, :)).^2);
    d = 1./d; % inverse so nearby channels have the largest weight
    d(badchans) = 0; % take out bad channels
    d = d./sum(d); % normalize so the total weight is one
    d = ascolumn(d);
    ISF.data(badchans(i), :) = sum(ISF.data .* repmat(d, 1, size(ISF.data, 2)), 1);
end
% Init the output structure
FIT = struct();
for i = 1:ISF.nbchan
    % Fit gaussians to the spectrum
    this_fit = fitisfspect('gaussianfit', ISF.freqs, ascolumn(ISF.data(i, :)));
    % Store spectral data and infraslow modulation peak
    FIT(i).freq = this_fit.freq;
    FIT(i).pow = this_fit.pow;
    FIT(i).fitpow = this_fit.fitpow;
    FIT(i).mu = this_fit.mu;
    FIT(i).peak = this_fit.peak;
    FIT(i).lower = this_fit.lower;
    FIT(i).upper = this_fit.upper;
    FIT(i).bandwidth = this_fit.bandwidth;
    FIT(i).coeffvals = this_fit.coeffvals;
    FIT(i).span = this_fit.span;
    FIT(i).offset = this_fit.offset;
end
% Extract the features used in group level analysis
ISF.features = css_extractfeatures(FIT, ISF);

end