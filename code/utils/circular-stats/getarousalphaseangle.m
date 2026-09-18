function [T] = getarousalphaseangle(SIG, const, filtcfg)
T = [];
cntbouts = struct();
cntbouts.total = 0;
cntbouts.tooshort = 0;
cntbouts.noarousal = 0;
cntbouts.processed = 0;
cntbouts.invalid = 0;
cntbouts.nplots = 0;
cntbouts.prearousal = 0;
% -------------------------------------------------------------------------
% Loop over files and extract bouts
for i = 1:length(SIG)
    % Get indices of the selected arousal bouts
    bouts = find(strcmpi({SIG(i).event.type}, 'boundary'));
    bouts = ceil([1, [SIG(i).event(bouts).latency], SIG(i).pnts]);
    % Normalization factor
    normfact = median(SIG(i).data, 2, 'omitnan');
    % Loop over bouts
    for b = 1:length(bouts)-1
        % Up the counter
        cntbouts.total = cntbouts.total+1;
        % ... skip if the bout is shorter than the filter order
        switch const.append_type
            case 'none'
                boutlength = (bouts(b+1)-bouts(b))+const.crop_aro(2);
            otherwise
                boutlength = (bouts(b+1)-bouts(b))+filtcfg.ford_mult*filtcfg.order;
        end
        % Check if bout is long enough
        if boutlength < 2*filtcfg.order+1
            cntbouts.tooshort = cntbouts.tooshort+1;
            continue
        end
        % Get the index of the selected event
        idx_evt = find([SIG(i).event.latency] > bouts(b+1)-30 * SIG(i).srate & [SIG(i).event.latency] < bouts(b+1)-29.8 * SIG(i).srate);
        idx_evt = idx_evt(contains({SIG(i).event(idx_evt).type}, 'arousal'));
        if isempty(idx_evt)
            cntbouts.noarousal = cntbouts.noarousal+1;
            continue
        end
        if length(idx_evt) ~= 1
            error('Found more than one selected arousal in bout ''%i''! In file ''%s''.', b, SIG(i).setname)
        end
        % Create filtered signal
        sig = pop_select(SIG(i), 'point', [bouts(b), bouts(b+1)]+const.crop_aro);
        sig.data = 10.*log10(sig.data./repmat(normfact, 1, sig.pnts));
        sig.times = sig.times./1000;
        % Exclude arousal events that have another arousal within 10 seconds of it
        idx_other = find(strcmpi({sig.event.type}, 'arousal') | strcmpi({sig.event.type}, 'arousalemg'));
        if ~isempty(idx_other)
            fin_aro = max([sig.event(idx_other).latency]+[sig.event(idx_other).duration]);
            if (sig.pnts - fin_aro) < (10*sig.srate-(abs(const.crop_aro(2))-30*sig.srate))
                cntbouts.prearousal = cntbouts.prearousal+1;
                continue
            end
        end
        % Up the counter
        cntbouts.processed = cntbouts.processed+1;
        % Calculate Hilbert on filtered signal (and predict last 5 seconds)
        fsig = isffilterbout(sig, filtcfg);
        [fsig, cmplx] = predict_isf(fsig, const);
        % Pivot the matrix
        instamp = abs(cmplx)'; % instantaneous amplitude
        phaseang = angle(cmplx)'; % instantaneous phase
        idx_aroonset = size(phaseang, 1);
        % Create new row for output table
        kv = filename2struct(SIG(i).setname);
        tmp = table();
        tmp.amp_d0 = ascolumn(instamp(idx_aroonset, :));
        tmp.phase_d0 = ascolumn(phaseang(idx_aroonset, :));
        tmp.phase_d1 = ascolumn(phaseang(idx_aroonset-1*fsig.srate, :));
        tmp.phase_d2 = ascolumn(phaseang(idx_aroonset-2*fsig.srate, :));
        tmp.phase_d3 = ascolumn(phaseang(idx_aroonset-3*fsig.srate, :));
        tmp.phase_d4 = ascolumn(phaseang(idx_aroonset-4*fsig.srate, :));
        tmp.phase_d8 = ascolumn(phaseang(idx_aroonset-8*fsig.srate, :));
        tmp.phase_d12 = ascolumn(phaseang(idx_aroonset-12*fsig.srate, :));
        tmp.phase_d16 = ascolumn(phaseang(idx_aroonset-16*fsig.srate, :));
        tmp.channel = {SIG(i).chanlocs.labels}';
        tmp.id = repmat({getuuid()}, size(tmp.phase_d0, 1), 1);
        tmp.sub = repmat({kv.sub}, size(tmp.phase_d0, 1), 1);
        tmp.ses = repmat({kv.ses}, size(tmp.phase_d0, 1), 1);
        tmp.bout = repmat(b, size(tmp.phase_d0, 1), 1);
        tmp.aro_type = repmat({SIG(i).event(idx_evt).type}, size(tmp.phase_d0, 1), 1);
        tmp.is_awakening = repmat(categorical(SIG(i).event(idx_evt).is_awakening), size(tmp.phase_d0, 1), 1);
        tmp.stage = repmat(SIG(i).event(idx_evt).stage, size(tmp.phase_d0, 1), 1);
        tmp.next_stage = repmat(SIG(i).event(idx_evt).next_stage, size(tmp.phase_d0, 1), 1);
        tmp.duration = repmat(SIG(i).event(idx_evt).duration, size(tmp.phase_d0, 1), 1);
        tmp.origlatency = repmat(SIG(i).event(idx_evt).origlatency, size(tmp.phase_d0, 1), 1);
        % Save data trace
        idx_this_aro_raw = (idx_aroonset-60*sig.srate+1:idx_aroonset)-(abs(const.crop_aro(2))-30*sig.srate);
        rawdata = sig.data(:, idx_this_aro_raw);
        tmp.rawdata = arrayfun(@(ridx) rawdata(ridx, :), 1:size(rawdata, 1), 'UniformOutput', false)';
        tmp.smtdata = arrayfun(@(ridx) asrow(smooth(rawdata(ridx, :), 100)), 1:size(rawdata, 1), 'UniformOutput', false)';
        isfdata = fsig.data(:, idx_aroonset-60*fsig.srate+1:idx_aroonset);
        tmp.isfdata = arrayfun(@(ridx) isfdata(ridx, :), 1:size(isfdata, 1), 'UniformOutput', false)';
        % Probability of other arousals
        idx_other = find(strcmpi({sig.event.type}, 'arousal') | strcmpi({sig.event.type}, 'arousalemg'));
        pr_aro = zeros(1, sig.pnts);
        idx_other = [[sig.event(idx_other).latency];[sig.event(idx_other).latency]+[sig.event(idx_other).duration]];
        if ~isempty(idx_other)
            idx_other = round(idx_other);
            idx_other = arrayfun(@(s,e) s:e, idx_other(1,:), idx_other(2,:), 'UniformOutput', false);
            idx_other = horzcat(idx_other{:});
            pr_aro(idx_other) = 1;
        end
        pr_aro = {pr_aro(idx_this_aro_raw)};
        tmp.pr_aro = repmat(pr_aro, size(tmp, 1), 1);
        % Save peaks and their locations
        [tmp.peak_amps, tmp.peak_idx] = zerocrosspeakfind(isfdata);
        % store the row in the output table
        if isempty(T)
            T = tmp;
        else
            T = [T; tmp]; %#ok<AGROW>
        end
    end
end
end