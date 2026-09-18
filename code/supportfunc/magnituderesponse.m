function [f, db, filtorder] = magnituderesponse(wintype, freqs, transbw, srate, filtorder)

% Filter order
if nargin < 5
    filtorder = pop_firwsord(lower(wintype), srate, transbw);
end

% Parse filter settings
c = [{filtorder} {sort(freqs / (srate / 2))}]; % Sorting and normalization
c = [c {windows(lower(wintype), filtorder + 1)'}];

% Get filter coefficients
b = firws(c{:});

nfft = 2^fix(log2(length(b)));
if nfft < 512
    nfft = 512;
end

% Frequency response
f = (0:1 / nfft:1) * srate / 2;
z = fft(b, nfft * 2);
z = z(1:fix(length(z) / 2) + 1);

% Magnitude response
db = abs(z);
db(db < eps^(2 / 3)) = eps^(2 / 3); % Log of zero warning
db = 20 * log10(db);

end