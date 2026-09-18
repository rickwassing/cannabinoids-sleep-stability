function [predProb, predCI] = calc_risk_ratio(ARO, mdl)
% Predictions
predTbl = ARO(1:2,:);
predTbl.cond(1) = {'placebo'};
predTbl.cond(2) = {'etc120'};
predTbl.duration = [0; 0];
predTbl.duration_z = [0; 0];

[predProb, predCI] = predict(mdl, predTbl, 'Conditional', false);

end
