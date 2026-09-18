function [f, phi, filtorder] = phaseresponse(wintype, freqs, transbw, srate, causal, filtorder)

% Filter order
if nargin < 6
    filtorder = pop_firwsord(lower(wintype), srate, transbw);
end

% Parse filter settings
c = [{filtorder} {sort(freqs / (srate / 2))}]; % Sorting and normalization
c = [c {windows(lower(wintype), filtorder + 1)'}];

% Get filter coefficients
b =  firws(c{:});

nfft = 2^fix(log2(length(b)));
if nfft < 512
    nfft = 512;
end

% Frequency response
f = (0:1 / nfft:1) * srate / 2;
m = fix((length(b) - 1) / 2); % Filter order
z = fft(b, nfft * 2);
z = z(1:fix(length(z) / 2) + 1);

% Phase response
z(abs(z) < eps^(2 / 3)) = NaN; % Phase is undefined for magnitude zero
phi = angle(z);
if causal
    phi = unwrap(phi);
else
    delay = -mod((0:1 / nfft:1) * m * pi + pi, 2 * pi) + pi; % Zero-phase
    phi = phi - delay;
    phi = phi + 2 * pi * (phi <= -pi + eps ^ (1/3)); % Unwrap
end

end