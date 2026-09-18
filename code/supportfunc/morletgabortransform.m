function [m_Transform, v_freq] = morletgabortransform(eegsignal, srate, freqlo, freqhi, freqres)

% This function calculates the wavelet transform of a signal of interest
% using a Gabor-Morlet kernel.
%
% Sintaxis:
%           f_MyGT(Y,Fr,minf,maxf,resf,b_2Draw)
% Inputs:
%           eegsignal: Signal to transform.
%           srate: Sample frequency
%           freqlo: min transform frequency. Default 10
%           freqhi: max transform frequency. Default 15
%           freqres: frequency resolution. Default 0.1
%
% Outputs:
%           m_Transform: Matrix containing the complex wavelet transform of
%           the signal between the designed frequencies. In order to get
%           the amplitude information use abs(m_Transform) and to get the
%           phase information use angle(m_Transform)
%
%           v_freq: vector containing the frequency information in
%           m_Transform.
%
% See also f_FFTPowAsLab, f_GeneralDynamics, f_PhaseAnalysis, f_SpindleDetection
%
% Alejandro Osorio-Forero 2019

% This is to put the signal as a column vector instead of a row vector
if size(eegsignal,1) < size(eegsignal,2)
    eegsignal = eegsignal';
end

cycles = 4;

v_freq = freqlo : freqres : freqhi;
s_Len = length(eegsignal);
s_HalfLen = floor(s_Len / 2) + 1;

v_YFFT = fft(eegsignal, numel(eegsignal));

v_WAxis = (2.* pi./ s_Len).* (0:(s_Len - 1)); % Vector with length of signal normalized to 2Pi
v_WAxis = v_WAxis.* srate;
v_WAxisHalf = v_WAxis(1:s_HalfLen); % first half of the signal, until Pi*sampling rate
m_Transform = zeros(numel(v_freq), length(eegsignal)); % Initial matrix

for iter = 1:length(v_freq)

    s_ActFrq = v_freq(iter); % Current Frequency

    dtseg = cycles * (1 /s_ActFrq); % Size of the gaussian window (pretty much the standard deviation of the gaussian)

    v_WinFFT = zeros(s_Len, 1);

    % Transform with normalizations
    v_WinFFT(1:s_HalfLen) = 2*exp(-0.5.* realpow(v_WAxisHalf - (2.* pi.* s_ActFrq), 2).*realpow(dtseg,2));
    v_WinFFT = v_WinFFT.* sqrt(s_Len)./ norm(v_WinFFT, 2); % Normalization L2
    
    % Normalization for unitary energy
    m_Transform(iter,:) = ifft(v_YFFT.* v_WinFFT)./ sqrt(dtseg);

end

end