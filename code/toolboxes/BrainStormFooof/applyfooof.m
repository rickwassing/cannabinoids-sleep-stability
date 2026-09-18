function fooof = applyfooof(Pow, Freq, varargin)

% 'Pow' must be a 2-D matrix of non-log transformed power <channels by frequecies>
% 'Freq' is a vector containing the frequencies

% Must be double and freq a row vector
Pow = double(Pow);
Freq = double(asrow(Freq));

% Define options
opt.freq_range = [0.1, 30];
opt.aperiodic_mode = 'fixed'; % 'fixed', 'knee', or 'floor'
opt.max_peaks = 3;
opt.peak_threshold = 2;
opt.min_peak_height = 3 / 10;
opt.peak_width_limits = [1.5, 12];
opt.proximity_threshold = 2;
opt.border_threshold = 1;
opt.peak_type = 'gaussian'; % 'gaussian', 'cauchy
opt.guess_weight = 'weak'; % 'none', 'weak', or 'strong'
opt.thresh_after = true;
opt.return_spectrum = 0;
opt.sort_type  = 'param'; % must be 'param'
opt.sort_param = 'frequency'; % 'frequency', 'amplitude' or 'std'

for i = 1:2:length(varargin)
    switch lower(varargin{i})
        case 'freq_range'
            opt.freq_range = varargin{i+1};
        case 'aperiodic_mode'
            opt.aperiodic_mode = varargin{i+1};
        case 'max_peaks'
            opt.max_peaks = varargin{i+1};
        case 'peak_threshold'
            opt.peak_threshold = varargin{i+1};
        case 'min_peak_height'
            opt.min_peak_height = varargin{i+1};
        case 'peak_width_limits'
            opt.peak_width_limits = varargin{i+1};
        case 'proximity_threshold'
            opt.proximity_threshold = varargin{i+1};
        case 'border_threshold'
            opt.border_threshold = varargin{i+1};
        case 'peak_type'
            opt.peak_type = varargin{i+1};
        case 'guess_weight'
            opt.guess_weight = varargin{i+1};
        case 'thresh_after'
            opt.thresh_after = varargin{i+1};
        case 'return_spectrum'
            opt.return_spectrum = varargin{i+1};
        case 'sort_type'
            opt.sort_type = varargin{i+1};
        case 'sort_param'
            opt.sort_param = varargin{i+1};
    end
end

% Boolean indicating whether we have access to optimization tools
hasOptimTools = false;

% Number of channels
nChan = size(Pow, 1);

% Find all frequency values within user limits
fMask = Freq >= opt.freq_range(1) & Freq <= opt.freq_range(2);
Freq = Freq(fMask);

% Log-transform power-spectrum
Pow = log10(squeeze(Pow(:, fMask))); % extract log spectra

% Initalize FOOOF structs
fooof(nChan) = struct(...
    'aperiodic_params', [],...
    'peak_params',      [],...
    'peak_types',       '',...
    'ap_fit',           [],...
    'fooofed_spectrum', [],...
    'peak_fit',         [],...
    'error',            [],...
    'r_squared',        []);

% Iterate across channels
for chan = 1:nChan
    % Fit aperiodic
    aperiodic_pars = robust_ap_fit(Freq, Pow(chan, :), opt.aperiodic_mode);
    % Remove aperiodic
    flat_spec = flatten_spectrum(Freq, Pow(chan, :), aperiodic_pars, opt.aperiodic_mode);
    % Fit peaks
    [peak_pars, peak_function] = fit_peaks(Freq, flat_spec, opt.max_peaks, opt.peak_threshold, opt.min_peak_height, ...
        opt.peak_width_limits/2, opt.proximity_threshold, opt.border_threshold, opt.peak_type, opt.guess_weight, hasOptimTools);
    if opt.thresh_after && ~hasOptimTools  % Check thresholding requirements are met for unbounded optimization
        peak_pars(peak_pars(:,2) < opt.min_peak_height,:)     = []; % remove peaks shorter than limit
        peak_pars(peak_pars(:,3) < opt.peak_width_limits(1)/2,:)  = []; % remove peaks narrower than limit
        peak_pars(peak_pars(:,3) > opt.peak_width_limits(2)/2,:)  = []; % remove peaks broader than limit
        peak_pars = drop_peak_cf(peak_pars, opt.border_threshold, opt.freq_range); % remove peaks outside frequency limits
        peak_pars(peak_pars(:,1) < 0,:) = []; % remove peaks with a centre frequency less than zero (bypass drop_peak_cf)
        peak_pars = drop_peak_overlap(peak_pars, opt.proximity_threshold); % remove smallest of two peaks fit too closely
    end
    % Refit aperiodic
    aperiodic = Pow(chan, :);
    for peak = 1:size(peak_pars,1)
        aperiodic = aperiodic - peak_function(Freq,peak_pars(peak,1), peak_pars(peak,2), peak_pars(peak,3));
    end
    aperiodic_pars = simple_ap_fit(Freq, aperiodic, opt.aperiodic_mode);
    % Generate model fit
    ap_fit = gen_aperiodic(Freq, aperiodic_pars, opt.aperiodic_mode);
    model_fit = ap_fit;
    for peak = 1:size(peak_pars,1)
        model_fit = model_fit + peak_function(Freq,peak_pars(peak,1),...
            peak_pars(peak,2),peak_pars(peak,3));
    end
    % Calculate model error
    MSE = sum((Pow(chan, :) - model_fit).^2)/length(model_fit);
    rsq_tmp = corrcoef(Pow(chan, :),model_fit).^2;
    % Return FOOOF results
    aperiodic_pars(2) = abs(aperiodic_pars(2));
    fooof(chan).aperiodic_params = aperiodic_pars;
    fooof(chan).peak_params      = peak_pars;
    fooof(chan).peak_types       = func2str(peak_function);
    fooof(chan).ap_fit           = 10.^ap_fit;
    fooof(chan).fooofed_spectrum = 10.^model_fit;
    fooof(chan).peak_fit         = 10.^(model_fit-ap_fit);
    fooof(chan).error            = MSE;
    fooof(chan).r_squared        = rsq_tmp(2);
    if opt.return_spectrum
        fooof(chan).power_spectrum = Pow(chan, :);
    end
    
end

end


% SUPPORT FUNCS

% ===== GENERATE APERIODIC =====
function ap_vals = gen_aperiodic(freqs,aperiodic_params,aperiodic_mode)
%       Generate aperiodic values, from parameter definition.
%
%       Parameters
%       ----------
%       freqs : 1xn array
%       	Frequency vector to create aperiodic component for.
%       aperiodic_params : 1x3 array
%           Parameters that define the aperiodic component.
%       aperiodic_mode : {'fixed', 'knee'}
%           Defines absence or presence of knee in aperiodic component.
%
%       Returns
%       -------
%       ap_vals : 1d array
%           Generated aperiodic values.
switch aperiodic_mode
    case 'fixed'  % no knee
        ap_vals = expo_nk_function(freqs,aperiodic_params);
    case 'knee'
        ap_vals = expo_function(freqs,aperiodic_params);
    case 'floor'
        ap_vals = expo_fl_function(freqs,aperiodic_params);
end
end

% ===== CORE MODELS =====
function ys = gaussian(freqs, mu, hgt, sigma)
%       Gaussian function to use for fitting.
%
%       Parameters
%       ----------
%       freqs : 1xn array
%           Frequency vector to create gaussian fit for.
%       mu, hgt, sigma : doubles
%           Parameters that define gaussian function (centre frequency,
%           height, and standard deviation).
%
%       Returns
%       -------
%       ys :    1xn array
%       Output values for gaussian function.
ys = hgt*exp(-(((freqs-mu)./sigma).^2) /2);
end

function ys = cauchy(freqs, ctr, hgt, gam)
%       Cauchy function to use for fitting.
%
%       Parameters
%       ----------
%       freqs : 1xn array
%           Frequency vector to create cauchy fit for.
%       ctr, hgt, gam : doubles
%           Parameters that define cauchy function (centre frequency,
%           height, and "standard deviation" [gamma]).
%
%       Returns
%       -------
%       ys :    1xn array
%       Output values for cauchy function.
ys = hgt./(1+((freqs-ctr)/gam).^2);
end

function ys = expo_function(freqs,params)
%       Exponential function to use for fitting 1/f, with a 'knee' (maximum at low frequencies).
%
%       Parameters
%       ----------
%       freqs : 1xn array
%           Input x-axis values.
%       params : 1x3 array (offset, knee, exp)
%           Parameters (offset, knee, exp) that define Lorentzian function:
%           y = 10^offset * (1/(knee + x^exp))
%
%       Returns
%       -------
%       ys :    1xn array
%           Output values for exponential function.
ys = params(1) - log10(abs(params(2)) +freqs.^params(3));
end

function ys = expo_nk_function(freqs, params)
%       Exponential function to use for fitting 1/f, without a 'knee'.
%
%       Parameters
%       ----------
%       freqs : 1xn array
%           Input x-axis values.
%       params : 1x2 array (offset, exp)
%           Parameters (offset, exp) that define Lorentzian function:
%           y = 10^offset * (1/(x^exp))
%
%       Returns
%       -------
%       ys :    1xn array
%           Output values for exponential (no-knee) function.
ys = params(1) - log10(freqs.^params(2));
end

function ys = expo_fl_function(freqs, params)
ys = log10(f.^(params(1)) * 10^(params(2)) + params(3));
end

% ===== FITTING ALGORITHM =====
function aperiodic_params = simple_ap_fit(freqs, power_spectrum, aperiodic_mode)
%       Fit the aperiodic component of the power spectrum.
%
%       Parameters
%       ----------
%       freqs : 1xn array
%           Frequency values for the power spectrum, in linear scale.
%       power_spectrum : 1xn array
%           Power values, in log10 scale.
%       aperiodic_mode : {'fixed','knee'}
%           Defines absence or presence of knee in aperiodic component.
%
%       Returns
%       -------
%       aperiodic_params : 1xn array
%           Parameter estimates for aperiodic fit.

%       Set guess params for lorentzian aperiodic fit, guess params set at init
options = optimset('Display', 'off', 'TolX', 1e-4, 'TolFun', 1e-6, ...
    'MaxFunEvals', 5000, 'MaxIter', 5000);

switch (aperiodic_mode)
    case 'fixed'  % no knee
        exp_guess = -(power_spectrum(end)-power_spectrum(1))./log10(freqs(end)./freqs(1));
        guess_vec = [power_spectrum(1), exp_guess];
        aperiodic_params = fminsearch(@error_expo_nk_function, guess_vec, options, freqs, power_spectrum);
    case 'knee'
        exp_guess = -(power_spectrum(end)-power_spectrum(1))./log10(freqs(end)./freqs(1));
        guess_vec = [power_spectrum(1),0, exp_guess];
        aperiodic_params = fminsearch(@error_expo_function, guess_vec, options, freqs, power_spectrum);
end

end

function aperiodic_params = robust_ap_fit(freqs, power_spectrum, aperiodic_mode)
%       Fit the aperiodic component of the power spectrum robustly, ignoring outliers.
%
%       Parameters
%       ----------
%       freqs : 1xn array
%           Frequency values for the power spectrum, in linear scale.
%       power_spectrum : 1xn array
%           Power values, in log10 scale.
%       aperiodic_mode : {'fixed','knee'}
%           Defines absence or presence of knee in aperiodic component.
%
%       Returns
%       -------
%       aperiodic_params : 1xn array
%           Parameter estimates for aperiodic fit.

% Do a quick, initial aperiodic fit
popt = simple_ap_fit(freqs, power_spectrum, aperiodic_mode);
initial_fit = gen_aperiodic(freqs, popt, aperiodic_mode);

% Flatten power_spectrum based on initial aperiodic fit
flatspec = power_spectrum - initial_fit;

% Flatten outliers - any points that drop below 0
flatspec(flatspec(:) < 0) = 0;

% Use percential threshold, in terms of # of points, to extract and re-fit
perc_thresh = prctile(flatspec, 0.025);
perc_mask = flatspec <= perc_thresh;
freqs_ignore = freqs(perc_mask);
spectrum_ignore = power_spectrum(perc_mask);

% Second aperiodic fit - using results of first fit as guess parameters
options = optimset('Display', 'off', 'TolX', 1e-4, 'TolFun', 1e-6, ...
    'MaxFunEvals', 5000, 'MaxIter', 5000);
guess_vec = popt;

switch (aperiodic_mode)
    case 'fixed'  % no knee
        aperiodic_params = fminsearch(@error_expo_nk_function, guess_vec, options, freqs_ignore, spectrum_ignore);
    case 'knee'
        aperiodic_params = fminsearch(@error_expo_function, guess_vec, options, freqs_ignore, spectrum_ignore);
end
end

function spectrum_flat = flatten_spectrum(freqs, power_spectrum, robust_aperiodic_params, aperiodic_mode)
%       Flatten the power spectrum by removing the aperiodic component.
%
%       Parameters
%       ----------
%       freqs : 1xn array
%           Frequency values for the power spectrum, in linear scale.
%       power_spectrum : 1xn array
%           Power values, in log10 scale.
%       robust_aperiodic_params : 1x2 or 1x3 array (see aperiodic_mode)
%           Parameter estimates for aperiodic fit.
%       aperiodic_mode : 1 or 2
%           Defines absence or presence of knee in aperiodic component.
%
%       Returns
%       -------
%       spectrum_flat : 1xn array
%           Flattened (aperiodic removed) power spectrum.

spectrum_flat = power_spectrum - gen_aperiodic(freqs,robust_aperiodic_params,aperiodic_mode);

end

function [model_params,peak_function] = fit_peaks(freqs, flat_iter, max_n_peaks, peak_threshold, min_peak_height, gauss_std_limits, proxThresh, bordThresh, peakType, guess_weight,hOT)
%       Iteratively fit peaks to flattened spectrum.
%
%       Parameters
%       ----------
%       freqs : 1xn array
%           Frequency values for the power spectrum, in linear scale.
%       flat_iter : 1xn array
%           Flattened (aperiodic removed) power spectrum.
%       max_n_peaks : double
%           Maximum number of gaussians to fit within the spectrum.
%       peak_threshold : double
%           Threshold (in standard deviations of noise floor) to detect a peak.
%       min_peak_height : double
%           Minimum height of a peak (in log10).
%       gauss_std_limits : 1x2 double
%           Limits to gaussian (cauchy) standard deviation (gamma) when detecting a peak.
%       proxThresh : double
%           Minimum distance between two peaks, in st. dev. (gamma) of peaks.
%       peakType : {'gaussian', 'cauchy', 'both'}
%           Which types of peaks are being fitted
%       guess_weight : {'none', 'weak', 'strong'}
%           Parameter to weigh initial estimates during optimization (None, Weak, or Strong)
%       hOT : 0 or 1
%           Defines whether to use constrained optimization, fmincon, or
%           basic simplex, fminsearch.
%
%       Returns
%       -------
%       gaussian_params : mx3 array, where m = No. of peaks.
%           Parameters that define the peak fit(s). Each row is a peak, as [mean, height, st. dev. (gamma)].
switch peakType
    case 'gaussian' % gaussian only
        peak_function = @gaussian; % Identify peaks as gaussian
        % Initialize matrix of guess parameters for gaussian fitting.
        guess_params = zeros(max_n_peaks, 3);
        % Save intact flat_spectrum
        flat_spec = flat_iter;
        % Find peak: Loop through, finding a candidate peak, and fitting with a guess gaussian.
        % Stopping procedure based on either the limit on # of peaks,
        % or the relative or absolute height thresholds.
        for guess = 1:max_n_peaks
            % Find candidate peak - the maximum point of the flattened spectrum.
            max_ind = find(flat_iter == max(flat_iter));
            max_height = flat_iter(max_ind);

            % Stop searching for peaks once max_height drops below height threshold.
            if max_height <= peak_threshold * std(flat_iter)
                break
            end

            % Set the guess parameters for gaussian fitting - mean and height.
            guess_freq = freqs(max_ind);
            guess_height = max_height;

            % Halt fitting process if candidate peak drops below minimum height.
            if guess_height <= min_peak_height
                break
            end

            % Data-driven first guess at standard deviation
            % Find half height index on each side of the center frequency.
            half_height = 0.5 * max_height;

            le_ind = sum(flat_iter(1:max_ind) <= half_height);
            ri_ind = length(flat_iter) - sum(flat_iter(max_ind:end) <= half_height)+1;

            % Keep bandwidth estimation from the shortest side.
            % We grab shortest to avoid estimating very large std from overalapping peaks.
            % Grab the shortest side, ignoring a side if the half max was not found.
            % Note: will fail if both le & ri ind's end up as None (probably shouldn't happen).
            short_side = min(abs([le_ind,ri_ind]-max_ind));
            % Estimate std from FWHM. Calculate FWHM, converting to Hz, get guess std from FWHM
            fwhm = short_side * 2 * (freqs(2)-freqs(1));
            guess_std = fwhm / (2 * sqrt(2 * log(2)));

            % Check that guess std isn't outside preset std limits; restrict if so.
            % Note: without this, curve_fitting fails if given guess > or < bounds.
            if guess_std < gauss_std_limits(1)
                guess_std = gauss_std_limits(1);
            end
            if guess_std > gauss_std_limits(2)
                guess_std = gauss_std_limits(2);
            end
            % Collect guess parameters.
            guess_params(guess,:) = [guess_freq, guess_height, guess_std];

            % Subtract best-guess gaussian.
            peak_gauss = gaussian(freqs, guess_freq, guess_height, guess_std);
            flat_iter = flat_iter - peak_gauss;

        end
        % Remove unused guesses
        guess_params(guess_params(:,1) == 0,:) = [];

        % Check peaks based on edges, and on overlap
        % Drop any that violate requirements.
        guess_params = drop_peak_cf(guess_params, bordThresh, [min(freqs) max(freqs)]);
        guess_params = drop_peak_overlap(guess_params, proxThresh);

        % If there are peak guesses, fit the peaks, and sort results.
        if ~isempty(guess_params)
            model_params = fit_peak_guess(guess_params, freqs, flat_spec, 1, guess_weight, gauss_std_limits,hOT);
        else
            model_params = zeros(1, 3);
        end

    case 'cauchy' % cauchy only
        peak_function = @cauchy; % Identify peaks as cauchy
        guess_params = zeros(max_n_peaks, 3);
        flat_spec = flat_iter;
        for guess = 1:max_n_peaks
            max_ind = find(flat_iter == max(flat_iter));
            max_height = flat_iter(max_ind);
            if max_height <= peak_threshold * std(flat_iter)
                break
            end
            guess_freq = freqs(max_ind);
            guess_height = max_height;
            if guess_height <= min_peak_height
                break
            end
            half_height = 0.5 * max_height;
            le_ind = sum(flat_iter(1:max_ind) <= half_height);
            ri_ind = length(flat_iter) - sum(flat_iter(max_ind:end) <= half_height);
            short_side = min(abs([le_ind,ri_ind]-max_ind));

            % Estimate gamma from FWHM. Calculate FWHM, converting to Hz, get guess gamma from FWHM
            fwhm = short_side * 2 * (freqs(2)-freqs(1));
            guess_gamma = fwhm/2;
            % Check that guess gamma isn't outside preset limits; restrict if so.
            % Note: without this, curve_fitting fails if given guess > or < bounds.
            if guess_gamma < gauss_std_limits(1)
                guess_gamma = gauss_std_limits(1);
            end
            if guess_gamma > gauss_std_limits(2)
                guess_gamma = gauss_std_limits(2);
            end

            % Collect guess parameters.
            guess_params(guess,:) = [guess_freq(1), guess_height, guess_gamma];

            % Subtract best-guess cauchy.
            peak_cauchy = cauchy(freqs, guess_freq(1), guess_height, guess_gamma);
            flat_iter = flat_iter - peak_cauchy;

        end
        guess_params(guess_params(:,1) == 0,:) = [];
        guess_params = drop_peak_cf(guess_params, bordThresh, [min(freqs) max(freqs)]);
        guess_params = drop_peak_overlap(guess_params, proxThresh);

        % If there are peak guesses, fit the peaks, and sort results.
        if ~isempty(guess_params)
            model_params = fit_peak_guess(guess_params, freqs, flat_spec, 2, guess_weight, gauss_std_limits,hOT);
        else
            model_params = zeros(1, 3);
        end
    case 'best' % best of both: model both fits and compare error, save best
        % Gaussian Fit
        guess_params = zeros(max_n_peaks, 3);
        flat_spec = flat_iter;
        for guess = 1:max_n_peaks
            max_ind = find(flat_iter == max(flat_iter));
            max_height = flat_iter(max_ind);
            if max_height <= peak_threshold * std(flat_iter)
                break
            end
            guess_freq = freqs(max_ind);
            guess_height = max_height;
            if guess_height <= min_peak_height
                break
            end
            half_height = 0.5 * max_height;
            le_ind = sum(flat_iter(1:max_ind) <= half_height);
            ri_ind = length(flat_iter) - sum(flat_iter(max_ind:end) <= half_height)+1;
            short_side = min(abs([le_ind,ri_ind]-max_ind));
            fwhm = short_side * 2 * (freqs(2)-freqs(1));
            guess_std = fwhm / (2 * sqrt(2 * log(2)));
            if guess_std < gauss_std_limits(1)
                guess_std = gauss_std_limits(1);
            end
            if guess_std > gauss_std_limits(2)
                guess_std = gauss_std_limits(2);
            end
            guess_params(guess,:) = [guess_freq, guess_height, guess_std];
            peak_gauss = gaussian(freqs, guess_freq, guess_height, guess_std);
            flat_iter = flat_iter - peak_gauss;
        end
        guess_params(guess_params(:,1) == 0,:) = [];
        guess_params = drop_peak_cf(guess_params, bordThresh, [min(freqs) max(freqs)]);
        guess_params = drop_peak_overlap(guess_params, proxThresh);
        if ~isempty(guess_params)
            gauss_params = fit_peak_guess(guess_params, freqs, flat_spec, 1, guess_weight, gauss_std_limits,hOT);
            flat_gauss = zeros(size(freqs));
            for peak = 1:size(gauss_params,1)
                flat_gauss =  flat_gauss + gaussian(freqs,gauss_params(peak,1),...
                    gauss_params(peak,2),gauss_params(peak,3));
            end
            error_gauss = sum((flat_gauss-flat_spec).^2);
        else
            gauss_params = zeros(1, 3); error_gauss = 1E10;
        end

        % Cauchy Fit
        guess_params = zeros(max_n_peaks, 3);
        flat_iter = flat_spec;
        for guess = 1:max_n_peaks
            max_ind = find(flat_iter == max(flat_iter));
            max_height = flat_iter(max_ind);
            if max_height <= peak_threshold * std(flat_iter)
                break
            end
            guess_freq = freqs(max_ind);
            guess_height = max_height;
            if guess_height <= min_peak_height
                break
            end
            half_height = 0.5 * max_height;
            le_ind = sum(flat_iter(1:max_ind) <= half_height);
            ri_ind = length(flat_iter) - sum(flat_iter(max_ind:end) <= half_height)+1;
            short_side = min(abs([le_ind,ri_ind]-max_ind));
            fwhm = short_side * 2 * (freqs(2)-freqs(1));
            guess_gamma = fwhm/2;
            if guess_gamma < gauss_std_limits(1)
                guess_gamma = gauss_std_limits(1);
            end
            if guess_gamma > gauss_std_limits(2)
                guess_gamma = gauss_std_limits(2);
            end
            guess_params(guess,:) = [guess_freq(1), guess_height, guess_gamma];
            peak_cauchy = cauchy(freqs, guess_freq(1), guess_height, guess_gamma);
            flat_iter = flat_iter - peak_cauchy;
        end
        guess_params(guess_params(:,1) == 0,:) = [];
        guess_params = drop_peak_cf(guess_params, bordThresh, [min(freqs) max(freqs)]);
        guess_params = drop_peak_overlap(guess_params, proxThresh);
        if ~isempty(guess_params)
            cauchy_params = fit_peak_guess(guess_params, freqs, flat_spec, 2, guess_weight, gauss_std_limits,hOT);
            flat_cauchy = zeros(size(freqs));
            for peak = 1:size(cauchy_params,1)
                flat_cauchy =  flat_cauchy + cauchy(freqs,cauchy_params(peak,1),...
                    cauchy_params(peak,2),cauchy_params(peak,3));
            end
            error_cauchy = sum((flat_cauchy-flat_spec).^2);
        else
            cauchy_params = zeros(1, 3); error_cauchy = 1E10;
        end
        % Save least-error model
        if min([error_gauss,error_cauchy]) == error_gauss
            model_params = gauss_params;
            peak_function = @gaussian;
        else
            model_params = cauchy_params;
            peak_function = @cauchy;
        end
end

end

function guess = drop_peak_cf(guess, bw_std_edge, freq_range)
%       Check whether to drop peaks based on center's proximity to the edge of the spectrum.
%
%       Parameters
%       ----------
%       guess : mx3 array, where m = No. of peaks.
%           Guess parameters for peak fits.
%
%       Returns
%       -------
%       guess : qx3 where q <= m No. of peaks.
%           Guess parameters for peak fits.

cf_params = guess(:,1)';
bw_params = guess(:,3)' * bw_std_edge;

% Check if peaks within drop threshold from the edge of the frequency range.

keep_peak = abs(cf_params-freq_range(1)) > bw_params & ...
    abs(cf_params-freq_range(2)) > bw_params;

% Drop peaks that fail the center frequency edge criterion
guess = guess(keep_peak,:);

end

function guess = drop_peak_overlap(guess, proxThresh)
%       Checks whether to drop gaussians based on amount of overlap.
%
%       Parameters
%       ----------
%       guess : mx3 array, where m = No. of peaks.
%           Guess parameters for peak fits.
%       proxThresh: double
%           Proximity threshold (in st. dev. or gamma) between two peaks.
%
%       Returns
%       -------
%       guess : qx3 where q <= m No. of peaks.
%           Guess parameters for peak fits.
%
%       Note
%       -----
%       For any gaussians with an overlap that crosses the threshold,
%       the lowest height guess guassian is dropped.

% Sort the peak guesses, so can check overlap of adjacent peaks
guess = sortrows(guess);

% Calculate standard deviation bounds for checking amount of overlap

bounds = [guess(:,1) - guess(:,3) * proxThresh, ...
    guess(:,1), guess(:,1) + guess(:,3) * proxThresh];

% Loop through peak bounds, comparing current bound to that of next peak
drop_inds =  [];

for ind = 1:size(bounds,1)-1

    b_0 = bounds(ind,:);
    b_1 = bounds(ind + 1,:);

    % Check if bound of current peak extends into next peak
    if b_0(2) > b_1(1)
        % If so, get the index of the gaussian with the lowest height (to drop)
        drop_inds = [drop_inds (ind - 1 + find(guess(ind:ind+1,2) == ...
            min(guess(ind,2),guess(ind+1,2))))];
    end
end
% Drop any peaks guesses that overlap too much, based on threshold.
guess(drop_inds,:) = [];
end

function peak_params = fit_peak_guess(guess, freqs, flat_spec, peak_type, guess_weight, std_limits, hOT)
%     Fits a group of peak guesses with a fit function.
%
%     Parameters
%     ----------
%       guess : mx3 array, where m = No. of peaks.
%           Guess parameters for peak fits.
%       freqs : 1xn array
%           Frequency values for the power spectrum, in linear scale.
%       flat_iter : 1xn array
%           Flattened (aperiodic removed) power spectrum.
%       peakType : {'gaussian', 'cauchy', 'best'}
%           Which types of peaks are being fitted.
%       guess_weight : 'none', 'weak', 'strong'
%           Parameter to weigh initial estimates during optimization.
%       std_limits: 1x2 array
%           Minimum and maximum standard deviations for distribution.
%       hOT : 0 or 1
%           Defines whether to use constrained optimization, fmincon, or
%           basic simplex, fminsearch.
%
%       Returns
%       -------
%       peak_params : mx3, where m =  No. of peaks.
%           Peak parameters post-optimization.


if hOT % Use OptimToolbox for fmincon
    options = optimset('Display', 'off', 'TolX', 1e-3, 'TolFun', 1e-5, ...
        'MaxFunEvals', 3000, 'MaxIter', 3000); % Tuned options
    lb = [guess(:,1)-guess(:,3)*2,zeros(size(guess(:,2))),ones(size(guess(:,3)))*std_limits(1)];
    ub = [guess(:,1)+guess(:,3)*2,inf(size(guess(:,2))),ones(size(guess(:,3)))*std_limits(2)];
    peak_params = fmincon(@error_model_constr,guess,[],[],[],[], ...
        lb,ub,[],options,freqs,flat_spec, peak_type);
else % Use basic simplex approach, fminsearch, with guess_weight
    options = optimset('Display', 'off', 'TolX', 1e-4, 'TolFun', 1e-5, ...
        'MaxFunEvals', 5000, 'MaxIter', 5000);
    peak_params = fminsearch(@error_model,...
        guess, options, freqs, flat_spec, peak_type, guess, guess_weight);
end
end

% ===== ERROR FUNCTIONS =====
function err = error_expo_nk_function(params,xs,ys)
ym = -log10(xs.^params(2)) + params(1);
err = sum((ys - ym).^2);
end

function err = error_expo_function(params,xs,ys)
ym = expo_function(xs,params);
err = sum((ys - ym).^2);
end

function err = error_model(params, xVals, yVals, peak_type, guess, guess_weight)
fitted_vals = 0;
weak = 1E2;
strong = 1E7;
for set = 1:size(params,1)
    switch (peak_type)
        case 'gaussian'
            fitted_vals = fitted_vals + gaussian(xVals, params(set,1), params(set,2), params(set,3));
        case 'cauchy'
            fitted_vals = fitted_vals + cauchy(xVals, params(set,1), params(set,2), params(set,3));
    end
end
switch guess_weight
    case 'none'
        err = sum((yVals - fitted_vals).^2);
    case 'weak' % Add small weight to deviations from guess m and amp
        err = sum((yVals - fitted_vals).^2) + ...
            weak*sum((params(:,1)-guess(:,1)).^2) + ...
            weak*sum((params(:,2)-guess(:,2)).^2);
    case 'strong' % Add large weight to deviations from guess m and amp
        err = sum((yVals - fitted_vals).^2) + ...
            strong*sum((params(:,1)-guess(:,1)).^2) + ...
            strong*sum((params(:,2)-guess(:,2)).^2);
end
end

function err = error_model_constr(params, xVals, yVals, peak_type)
fitted_vals = 0;
for set = 1:size(params,1)
    switch (peak_type)
        case 'gaussian'
            fitted_vals = fitted_vals + gaussian(xVals, params(set,1), params(set,2), params(set,3));
        case 'cauchy'
            fitted_vals = fitted_vals + cauchy(xVals, params(set,1), params(set,2), params(set,3));
    end
end
err = sum((yVals - fitted_vals).^2);
end


% ===== FOOOF STATS =====
function [ePeaks, eAper, eStats] = FOOOF_analysis(FOOOF_data, ChanNames, TF, max_peaks, sort_type, sort_param, sort_bands)
% Organize/extract peak components from FOOOF models
nChan = numel(ChanNames);
maxEnt = nChan * max_peaks;
switch sort_type
    case 'param'
        % Initialize output struct
        ePeaks = struct('channel', [], 'center_frequency', [],...
            'amplitude', [], 'std_dev', []);
        % Collect data from all peaks
        i = 0;
        for chan = 1:nChan
            if ~isempty(FOOOF_data(chan).peak_params)
                for p = 1:size(FOOOF_data(chan).peak_params,1)
                    i = i +1;
                    ePeaks(i).channel = ChanNames(chan);
                    ePeaks(i).center_frequency = FOOOF_data(chan).peak_params(p,1);
                    ePeaks(i).amplitude = FOOOF_data(chan).peak_params(p,2);
                    ePeaks(i).std_dev = FOOOF_data(chan).peak_params(p,3);
                end
            end
        end
        % Apply specified sort
        switch sort_param
            case 'frequency'
                [tmp,iSort] = sort([ePeaks.center_frequency]);
                ePeaks = ePeaks(iSort);
            case 'amplitude'
                [tmp,iSort] = sort([ePeaks.amplitude]);
                ePeaks = ePeaks(iSort(end:-1:1));
            case 'std'
                [tmp,iSort] = sort([ePeaks.std_dev]);
                ePeaks = ePeaks(iSort);
        end
    case 'band'
        % Initialize output struct
        ePeaks = struct('channel', [], 'center_frequency', [],...
            'amplitude', [], 'std_dev', [], 'band', []);
        % Generate bands from input
        bands = process_tf_bands('Eval', sort_bands);
        % Collect data from all peaks
        i = 0;
        for chan = 1:nChan
            if ~isempty(FOOOF_data(chan).peak_params)
                for p = 1:size(FOOOF_data(chan).peak_params,1)
                    i = i +1;
                    ePeaks(i).channel = ChanNames(chan);
                    ePeaks(i).center_frequency = FOOOF_data(chan).peak_params(p,1);
                    ePeaks(i).amplitude = FOOOF_data(chan).peak_params(p,2);
                    ePeaks(i).std_dev = FOOOF_data(chan).peak_params(p,3);
                    % Find name of frequency band from user definitions
                    bandRanges = cell2mat(bands(:,2));
                    iBand = find(ePeaks(i).center_frequency >= bandRanges(:,1) & ePeaks(i).center_frequency <= bandRanges(:,2));
                    if ~isempty(iBand)
                        ePeaks(i).band = bands{iBand,1};
                    else
                        ePeaks(i).band = 'None';
                    end
                end
            end
        end
end

% ===== EXTRACT APERIODIC =====
% Organize/extract aperiodic components from FOOOF models
hasKnee = length(FOOOF_data(1).aperiodic_params) - 2;
% Initialize output struct
eAper = struct('channel', [], 'offset', [], 'exponent', []);
for chan = 1:nChan
    eAper(chan).channel = ChanNames(chan);
    eAper(chan).offset = FOOOF_data(chan).aperiodic_params(1);
    if hasKnee % Legacy FOOOF alters order of parameters
        eAper(chan).exponent = FOOOF_data(chan).aperiodic_params(3);
        eAper(chan).knee_frequency = FOOOF_data(chan).aperiodic_params(2);
    else
        eAper(chan).exponent = FOOOF_data(chan).aperiodic_params(2);
    end
end

% ===== EXTRACT STAT =====
% Organize/extract stats from FOOOF models
% Initialize output struct
eStats = struct('channel', ChanNames);
for chan = 1:nChan
    eStats(chan).MSE = FOOOF_data(chan).error;
    eStats(chan).r_squared = FOOOF_data(chan).r_squared;
    spec = squeeze(log10(TF(chan,:,:)));
    fspec = squeeze(log10(FOOOF_data(chan).fooofed_spectrum))';
    eStats(chan).frequency_wise_error = abs(spec-fspec);
end
end
