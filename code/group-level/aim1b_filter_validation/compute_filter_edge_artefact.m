function [SIGA, SIGB, ANGA, ANGB, DANG, DAMP, sig_control, sig_test] = compute_filter_edge_artefact(SIGMA, filtcfg, append_type, ford_mult)
% -------------------------------------------------------------------------
% Create two copies of the sigma-timeseries. One which is 6 times the
% filter length, and another that is 4 times the filter length (cropped
% from the end). Then filter these timeseries, crop the long timeseries to
% the same size as the shorter one, and apply Hilbert transform to get the
% instantaneous phase estimate. The difference in the phase estimate is the
% impact of the edge artifact of the filter.
% -------------------------------------------------------------------------
% Init
H = struct();
ANGA = [];
ANGB = [];
SIGA = [];
SIGB = [];
DANG = [];
DAMP = [];
% -------------------------------------------------------------------------
% Loop over files and extract bouts
cbouts = struct();
cbouts.total = 0;
cbouts.tooshort = 0;
cbouts.analysed = 0;
trem = now(); %#ok<TNOW1>
for i = 1:length(SIGMA)
    % ---------------------------------------------------------------------
    % Get indices of the selected NREM bouts
    bouts = find(strcmpi({SIGMA(i).event.type}, 'boundary'));
    bouts = ceil([1, [SIGMA(i).event(bouts).latency], SIGMA(i).pnts]);
    % ---------------------------------------------------------------------
    % Loop over bouts...
    for b = 1:length(bouts)-1
        % ... continue to the next bout if is too short
        cbouts.total = cbouts.total+1;
        if bouts(b+1)-bouts(b) < 6*filtcfg.order
            cbouts.tooshort = cbouts.tooshort+1;
            continue
        end
        cbouts.analysed = cbouts.analysed+1;
        % -----------------------------------------------------------------
        % Create control signal
        sig_control = pop_select(SIGMA(i), 'point', [bouts(b), bouts(b+1)]);
        sig_control.times = sig_control.times./1000;
        sig_control.data = detrend(sig_control.data', 0, 'omitnan')'; % demean
        % -----------------------------------------------------------------
        % Create test signal
        sig_test = pop_select(SIGMA(i), 'point', [bouts(b), bouts(b+1)]+[0 -ford_mult*filtcfg.order]);
        sig_test.times = sig_test.times./1000;
        sig_test.data = detrend(sig_test.data', 0, 'omitnan')'; % demean
        % -----------------------------------------------------------------
        % Get prepend and append signals
        sig_test = signalappend(sig_test, append_type, ford_mult, filtcfg);
        sig_control = signalappend(sig_control, append_type, ford_mult, filtcfg);
        % -----------------------------------------------------------------
        % Prepend and append data to cover for filter edge artefact
        sig_test.data = [sig_test.prepend, sig_test.data, sig_test.append];
        sig_test.pnts = size(sig_test.data, 2);
        sig_test.xmin = 0;
        sig_test.xmax = sig_test.pnts/sig_test.srate;
        sig_test.times = linspace(sig_test.xmin, sig_test.xmax, sig_test.pnts);
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        sig_control.data = [sig_control.prepend, sig_control.data, sig_control.append];
        sig_control.pnts = size(sig_control.data, 2);
        sig_control.xmin = 0;
        sig_control.xmax = sig_control.pnts/sig_control.srate;
        sig_control.times = linspace(sig_control.xmin, sig_control.xmax, sig_control.pnts);
        % -----------------------------------------------------------------
        % Filter the data
        sig_control = pop_firws(sig_control, ...
            'fcutoff', filtcfg.cutoff+filtcfg.adj, ...
            'ftype', 'bandpass', ...
            'wtype', filtcfg.wintype, ...
            'warg', filtcfg.warg, ...
            'forder', filtcfg.order, ...
            'plotfresp', false, ...
            'minphase', 0);
        sig_test = pop_firws(sig_test, ...
            'fcutoff', filtcfg.cutoff+filtcfg.adj, ...
            'ftype', 'bandpass', ...
            'wtype', filtcfg.wintype, ...
            'warg', filtcfg.warg, ...
            'forder', filtcfg.order, ...
            'plotfresp', false, ...
            'minphase', 0);
        % -----------------------------------------------------------------
        % cut away the appended data
        sig_control.data = sig_control.data(:, length(sig_control.prepend)+1:end-length(sig_control.append));
        sig_control.pnts = size(sig_control.data, 2);
        sig_control.xmin = 0;
        sig_control.xmax = sig_control.pnts/sig_control.srate;
        sig_control.times = linspace(sig_control.xmin, sig_control.xmax, sig_control.pnts);
        % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
        sig_test.data = sig_test.data(:, length(sig_test.prepend)+1:end-length(sig_test.append));
        sig_test.pnts = size(sig_test.data, 2);
        sig_test.xmin = 0;
        sig_test.xmax = sig_test.pnts/sig_test.srate;
        sig_test.times = linspace(sig_test.xmin, sig_test.xmax, sig_test.pnts);
        % -----------------------------------------------------------------
        % Apply Hilbert to control signal
        H.control.sig = sig_control.data(:, 1:end-ford_mult*filtcfg.order)';
        H.control.sig = detrend(H.control.sig, 0);
        H.control.x = hilbert(H.control.sig);
        H.control.ang = angle(H.control.x);
        H.control.amp = abs(H.control.x);
        % Apply Hilbert to test signal
        H.test.sig = sig_test.data';
        H.test.sig = detrend(H.test.sig, 0);
        H.test.x = hilbert(H.test.sig);
        H.test.ang = angle(H.test.x);
        H.test.amp = abs(H.test.x);
        % Caclulate the difference in phase and amplitude
        dang = circ_dist(H.control.ang, H.test.ang);
        damp = abs(H.control.amp - H.test.amp);
        % Store estimated angle and difference in angle for this bout
        SIGA = [SIGA; H.control.sig(end-4*filtcfg.order+1:end, :)']; %#ok<AGROW>
        SIGB = [SIGB; H.test.sig(end-4*filtcfg.order+1:end, :)']; %#ok<AGROW>
        ANGA = [ANGA; H.control.ang(end-4*filtcfg.order+1:end, :)']; %#ok<AGROW>
        ANGB = [ANGB; H.test.ang(end-4*filtcfg.order+1:end, :)']; %#ok<AGROW>
        DANG = [DANG; dang(end-4*filtcfg.order+1:end, :)']; %#ok<AGROW>
        DAMP = [DAMP; damp(end-4*filtcfg.order+1:end, :)']; %#ok<AGROW>
    end
    trem = remainingTime(trem, length(SIGMA), true);
end

end
