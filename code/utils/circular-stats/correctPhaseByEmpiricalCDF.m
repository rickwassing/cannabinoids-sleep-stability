function theta_corr = correctPhaseByEmpiricalCDF(theta)
% correctPhaseByEmpiricalCDF
%
% Corrects circular phase angles for a non-uniform background phase
% distribution using a probability integral transform (PIT).
%
% INPUTS
%   theta     : vector of phase angles (radians, [-pi pi])
%   binEdges  : histogram bin edges spanning [-pi pi]
%   pdfVals   : PDF values per bin (from histcounts, 'Normalization','pdf')
%
% OUTPUT
%   theta_corr : phase angles transformed to be uniform under the null
%
% This transformation preserves circular geometry and allows standard
% circular non-uniformity tests to be applied relative to an empirical null.

% Ensure column vectors
theta = theta(:);
pdfVals = [0.168601227889208   0.164595880397698   0.158654764913127   0.154324463102062   0.152130761896888   0.152181056566155   0.152730375866730   0.154655992642140  0.156712767763782   0.160673011549596   0.165711706856006   0.168887307659352];
pdfVals = pdfVals(:);

Resolution = 30;
nbins = 360/Resolution;
binEdges = linspace(-pi, pi, nbins+1);

% Bin width (assumed constant)
binWidth = mean(diff(binEdges));

% Convert PDF to probability mass per bin
Pbin = pdfVals * binWidth;

% Empirical CDF defined at bin edges
CDF = [0; cumsum(Pbin)];

% Wrap phases to [-pi pi]
if any(theta < -pi) || any(theta > pi)
    error('theta must be in the range if [-pi, pi]');
end

% Probability integral transform
U = interp1(binEdges(:), CDF, theta, 'linear', 'extrap');

% Map back to circular domain
theta_corr = 2*pi*U - pi;
end