function Clusters = getpermpvalue(perms, fld, cond)

if isempty(cond)
    % Extract the permuation max-lengths

    XData = -90:0.1:-0.1;
    maxLengths = arrayfun(@(s) max([0, find(diff([0, s.(fld).pval < 0.05, 0]) == -1) - find(diff([0, s.(fld).pval < 0.05, 0]) == 1)]), perms, 'UniformOutput', false);
    maxLengths = cat(1, maxLengths{:});

    Clusters = struct();
    Clusters.idx = [find(diff([0, perms(1).(fld).pval < 0.05, 0]) == 1); find(diff([0, perms(1).(fld).pval < 0.05, 0]) == -1)];
    Clusters.ts = [...
        XData(Clusters.idx(1,:));
        XData(Clusters.idx(2,:))];
    Clusters.length = diff(Clusters.idx);
    for i = 1:length(Clusters.length)
        Clusters.perm_p(i) = sum(maxLengths >= Clusters.length(i))/length(maxLengths);
    end

    return
end

% Extract the permuation max-lengths

XData = -61.5:0.1:-1.5;
maxLengths = arrayfun(@(s) max([0, find(diff([0, s.(fld).(cond).pval < 0.05, 0]) == -1) - find(diff([0, s.(fld).(cond).pval < 0.05, 0]) == 1)]), perms, 'UniformOutput', false);
maxLengths = cat(1, maxLengths{:});

Clusters = struct();
Clusters.idx = [find(diff([0, perms(1).(fld).(cond).pval < 0.05, 0]) == 1); find(diff([0, perms(1).(fld).(cond).pval < 0.05, 0]) == -1)];
Clusters.ts = [...
    XData(Clusters.idx(1,:));
    XData(Clusters.idx(2,:))];
Clusters.length = diff(Clusters.idx);
for i = 1:length(Clusters.length)
    Clusters.perm_p(i) = sum(maxLengths >= Clusters.length(i))/length(maxLengths);
end

end