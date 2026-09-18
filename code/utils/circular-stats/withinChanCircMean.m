function MU = withinChanCircMean(T, fld)

if iscell(T.(fld))
    Data = cat(1, T.(fld){:});
else
    Data = cat(1, T.(fld)(:));
end
Evs = unique(T.id);
MU = nan(length(Evs), size(Data, 2));

for s = 1:length(Evs)
    idx_ev = strcmpi(T.id, Evs{s});
    T_ev = T(idx_ev, :);

    if iscell(T_ev.(fld))
        D_ev = cat(1, T_ev.(fld){:});
    else
        D_ev = cat(1, T_ev.(fld)(:));
    end

    MU(s, :) = circ_mean(D_ev);
end

end