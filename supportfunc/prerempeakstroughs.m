function [event] = prerempeakstroughs(h, cfg)

event = table();

for chan = 1:size(h.sigma, 2)
    [~, idx_rem_onset] = min(abs(h.times));
    fpksig = [0; abs(diff(h.ang(:, chan))) > cfg.pkheight];
    [~, idx_tr_all] = findpeaks(fpksig, 'MinPeakDistance', cfg.pkdist*cfg.srate);

    % Find the last previous trough and subsequent peak
    idx_sel = find(idx_tr_all <= idx_rem_onset - cfg.mindelay*cfg.srate, 1, 'last');
    if idx_sel+1 > length(idx_tr_all)
        idx_tr_rem = [idx_tr_all(idx_sel); min([idx_tr_all(idx_sel)+72.5*cfg.srate, size(h.sigma, 1)])];
    else
        idx_tr_rem = idx_tr_all([idx_sel, idx_sel+1]);
    end
    [~, idx_pk_rem] = max(h.sigma(idx_tr_rem(1):idx_tr_rem(2), chan));
    idx_pk_rem = idx_pk_rem+idx_tr_rem(1);
    idx_tr_rem = idx_tr_rem(1);

    % find previous troughs and peaks (not associated with REM transition)
    idx_tr_nrem = [idx_tr_all(idx_tr_all < idx_tr_rem); idx_tr_rem];
    idx_pk_nrem = [];
    for j = 1:length(idx_tr_nrem)-1
        [~, idx_pk_nrem(j)] = max(h.sigma(idx_tr_nrem(j):idx_tr_nrem(j+1), chan)); %#ok<AGROW>
        idx_pk_nrem(j) = idx_pk_nrem(j)+idx_tr_nrem(j); %#ok<AGROW>
    end
    idx_tr_nrem = idx_tr_nrem(1:end-1);
    amp_nrem = [];
    dur_nrem = [];
    for j = 1:length(idx_tr_nrem)
        amp_nrem = [amp_nrem; h.sigma(idx_pk_nrem(j), chan) - h.sigma(idx_tr_nrem(j), chan)]; %#ok<AGROW>
        dur_nrem = [dur_nrem; idx_pk_nrem(j) - idx_tr_nrem(j)]; %#ok<AGROW>
    end
    ang_rem = h.ang(idx_rem_onset, chan);
    amp_rem = h.sigma(idx_pk_rem, chan) - h.sigma(idx_tr_rem, chan);
    dur_rem = idx_pk_rem - idx_tr_rem;
    % Extract the ISF parameters of the prior NREM period
    tmp = table();
    tmp.latency = ascolumn(idx_tr_nrem);
    tmp.peak_latency = ascolumn(idx_pk_nrem);
    tmp.duration = ascolumn(dur_nrem);
    tmp.amplitude = ascolumn(amp_nrem);
    tmp.angle = nan(length(idx_tr_nrem), 1);
    tmp.idx_tr = repmat({nan}, length(idx_tr_nrem), 1);
    tmp.stage = repmat({'nrem'}, length(idx_tr_nrem), 1);
    tmp.channel = repmat(chan, length(idx_tr_nrem), 1);
    % Append
    event = [event; tmp]; %#ok<AGROW>
    % Extract the ISF parameters at the REM transition
    tmp = table();
    tmp.latency = ascolumn(idx_tr_rem);
    tmp.peak_latency = ascolumn(idx_pk_rem);
    tmp.duration = ascolumn(dur_rem);
    tmp.amplitude = ascolumn(amp_rem);
    tmp.angle = ang_rem;
    tmp.idx_tr = {idx_tr_all};
    tmp.stage = repmat({'rem'}, length(idx_tr_rem), 1);
    tmp.channel = repmat(chan, length(idx_tr_rem), 1);
    % Append
    event = [event; tmp]; %#ok<AGROW>
end

end