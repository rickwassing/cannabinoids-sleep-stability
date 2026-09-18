function [MU, EV] = withinSubMean(T, fld, smoothfactor)

if nargin < 3
    smoothfactor = 1;
end

if iscell(T.(fld))
    Data = cat(1, T.(fld){:});
else
    Data = cat(1, T.(fld)(:));
end
Subs = unique(T.sub);
MU = nan(length(Subs), size(Data, 2));
EV = [];

for s = 1:length(Subs)
    idx_sub = strcmpi(T.sub, Subs{s});
    Tsub = T(idx_sub, :);

    if iscell(Tsub.(fld))
        Dsub = cat(1, Tsub.(fld){:});
    else
        Dsub = cat(1, Tsub.(fld)(:));
    end

    if smoothfactor > 1
        for i = 1:size(Dsub, 1)
            Dsub(i, :) = smooth(Dsub(i, :), smoothfactor);
        end
    end

    ids = unique(Tsub.id);
    mu = nan(length(ids), size(Dsub, 2));
    for i = 1:length(ids)
        idx_ev = strcmpi(Tsub.id, ids{i});
        mu(i, :) = mean(Dsub(idx_ev, :), 1, 'omitnan');
    end
    MU(s, :) = mean(mu, 1, 'omitnan');
    EV = [EV; mu];
end

end