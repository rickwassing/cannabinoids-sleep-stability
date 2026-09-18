function normdata = withinConditionNorm(T, method)

Subs = unique(T.sub);
Conds = unique(T.ses);
normdata = [];
rt = now();
for i = 1:length(Subs)
    for j = 1:length(Conds)
        idx = find(strcmpi(T.sub, Subs{i}) & strcmpi(T.ses, Conds{j}));
        rawdata = T.rawdata(idx);
        switch lower(method)
            case 'db'
                tmp = cellfun(@(ts) 10.*log10(ts./median(ts)), rawdata, 'UniformOutput', false);
        end
        normdata = [normdata;tmp]; %#ok<AGROW> 
    end
    rt = remainingTime(rt, length(Subs));
end

end