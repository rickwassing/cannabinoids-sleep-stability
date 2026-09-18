function [avsigma, desmat] = build_rem_transition_features(H)
% -------------------------------------------------------------------------
% For each trace in 'H', extract the average raw sigma-power timeseries and
% the per-event peak amplitude/half-duration features used for modelling.
% -------------------------------------------------------------------------
avsigma = [];
desmat = [];
for i = 1:length(H)

    tmp = table();
    tmp.sub = {H(i).sub};
    tmp.cond = {H(i).cond};
    tmp.s = {mean(H(i).rawsigma(1:3601, :), 2)};
    avsigma = [avsigma; tmp]; %#ok<AGROW>

    events = H(i).event;

    uuids = ascolumn(sort(arrayfun(@(x) getuuid(), 1:50, 'UniformOutput', false)));
    [eventGroups, eventKeys] = findgroups(events.id); %#ok<ASGLU>

    tmp = table();
    tmp.latency = splitapply(@mean, events.latency+events.duration/2, eventGroups);
    tmp.duration = splitapply(@mean, events.duration/10, eventGroups);
    tmp.amplitude = splitapply(@mean, double(events.amplitude), eventGroups);
    tmp.stage = splitapply(@(x) x(1), events.stage, eventGroups);
    tmp.id = splitapply(@(x) x(1), events.id, eventGroups);
    tmp.sub = repmat({H(i).sub}, size(tmp, 1), 1);
    tmp.cond = repmat({H(i).cond}, size(tmp, 1), 1);
    tmp.remepdur = repmat(H(i).remepdur, size(tmp, 1), 1)./(60*10);
    tmp.remeplat = repmat(H(i).remeplat, size(tmp, 1), 1)./(60*10);
    tmp.id = uuids(findgroups(tmp.remepdur));

    desmat = [desmat; tmp]; %#ok<AGROW>
end

desmat.dur_z = zscore(desmat.duration);
desmat.amp_z = zscore(desmat.amplitude);
desmat.remdur_z = zscore(desmat.remepdur);
desmat.remlat_z = zscore(desmat.remeplat);

end
