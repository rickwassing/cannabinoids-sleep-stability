% NOTE: no current call sites found in code/ as of 2026-09-18 (Phase 3
% refactor). Kept in utils/misc/ rather than archived, since unclear
% current usage does not necessarily mean unused; flagged for your review.
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