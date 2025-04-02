function f = fitisfspect(fitMethod, freq, pow)

% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Fit gaussian curves, or use 'smoothdata'
switch fitMethod

    case 'gaussianfit'
        ngaussians = 2;
        fitpow = ascolumn(pow);
        [fitpow, offset] = gaussianfit(freq, fitpow, ngaussians);
        coeffvals = coeffvalues(fitpow);
        fitpow = fitpow(freq) + offset;
        span = nan;
        mu = coeffvals(2);
        [~, idx_mu] = min(abs(freq - mu));
        pks = fitpow(idx_mu);
        bandwidth = coeffvals(3);
        lolim = mu - bandwidth/2;
        hilim = mu + bandwidth/2;

    case 'smoothdata'
        span = 1;
        fitpow = smoothdata(pow, 'gaussian', span);
        pks = findpeaks(fitpow);
        while length(pks) > 1 && span < length(pow)
            span = span + 1;
            fitpow = smoothdata(pow, 'gaussian', span);
            pks = findpeaks(fitpow(fitpow > mean(fitpow)));
        end
        fitpow(freq < 0.0075) = 0;
        [pks, idx_max] = max(fitpow);
        mu = freq(idx_max);
        idx_edge = idx_max;
        while fitpow(idx_edge) > fitpow(idx_edge+1)
            idx_edge = idx_edge+1;
            if idx_edge == length(fitpow)
                break
            end
        end
        lolim = mu - (freq(idx_edge) - mu);
        hilim = freq(idx_edge);
        bandwidth = hilim - lolim;
        coeffvals = NaN;
        offset = 0;
        
    otherwise
        error('no such method.')
end

f = struct();
f.freq = freq;
f.pow = pow;
f.fitpow = fitpow;
f.mu = mu;
f.peak = pks;
f.lower = lolim;
f.upper = hilim;
f.bandwidth = bandwidth;
f.coeffvals = coeffvals;
f.span = span;
f.offset = offset;

end