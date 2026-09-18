% ARCHIVED 2026-09-18: no call sites found anywhere in code/ (SRP_REFACTOR_PLAN.md
% Phase 1 dead-code audit); kept for provenance, not called from main.m.
% -------------------------------------------------------------------------
function data = zscoreacrosschannels(data)
if size(data, 2) < size(data, 1)
    data = data'; % make sure the columns are channels
end
mu = data(:);
sd = std(mu);
mu = mean(mu);
data = (data-mu)./sd;

end