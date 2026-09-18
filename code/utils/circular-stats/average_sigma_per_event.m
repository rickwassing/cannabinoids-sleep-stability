% AVERAGE_SIGMA_PER_EVENT Average sigma power across all 18 channels for each arousal event
%
% INPUT:
%   T - table with columns: sub, cond, id, chan, Sigma, outcome (and optionally others)
%
% OUTPUT:
%   T_event - table with one row per event, with averaged Sigma across channels

function desmat_event = average_sigma_per_event(desmat)

% Ensure categorical variables
desmat.cond = categorical(desmat.cond);
desmat.is_awakening = categorical(desmat.is_awakening);
desmat.sub = categorical(desmat.sub);
desmat.id = categorical(desmat.id);

% Find unique events
[eventGroups, eventKeys] = findgroups(desmat.id);

% Average sigma across channels per event
y_avg = nan(length(eventKeys), 600);
pr_avg = nan(length(eventKeys), 600);
for t = 1:600
    if t == 1
        y_amp = splitapply(@mean, desmat.y_amp, eventGroups);
    end
    y_avg(:, t) = splitapply(@mean, desmat.(sprintf('y%i', t)), eventGroups);
    pr_avg(:, t) = splitapply(@mean, desmat.(sprintf('Pr%i', t)), eventGroups);
end

% Build one-row-per-event table
desmat_event = [...
    array2table(y_amp, 'VariableNames', "y_amp"), ...
    array2table(y_avg, 'VariableNames', strcat("y", string(1:600))), ...
    array2table(pr_avg, 'VariableNames', strcat("Pr", string(1:600)))];

% Recover subject and condition for each event
desmat_event.cond = splitapply(@(x) x(1), desmat.cond, eventGroups);
desmat_event.is_awakening = splitapply(@(x) x(1), desmat.is_awakening, eventGroups);
desmat_event.sub = splitapply(@(x) x(1), desmat.sub, eventGroups);
desmat_event.id = splitapply(@(x) x(1), desmat.id, eventGroups);

desmat_event = sortrows(desmat_event, {'sub', 'cond', 'is_awakening'});


end