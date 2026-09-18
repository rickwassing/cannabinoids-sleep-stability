function bouts = getnrembouts(HYP, srate, cutoff)
% Output variable bouts are in seconds, cutoff must also be in seconds
% -------------------------------------------------------------------------
% Initialize
bouts = [];
start = false;
cnt = 0; % Counter to keep track of how many non-stage 2 epochs occurred in a bout
% -------------------------------------------------------------------------
% For each epoch in the sleep staging
for i = 1:size(HYP, 1)
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Start a new bout at the first occurance of Stage 2 NREM sleep
    if ~start && HYP.stage_num(i) == -2
        start = true; % A new Stage 2 NREM bout has started
        thisbout = round(HYP.latency(i)); % Store the onset
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % If the current stage is 2, then reset the counter
    if HYP.stage_num(i) == -2
        cnt = 0;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % If the current stage in a bout is not 2, then increment the counter
    if start && HYP.stage_num(i) ~= -2
        cnt = cnt+1;
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % If there were more than 1 epochs of non-stage 2 then end this bout
    if start && cnt > 1
        start = false;
        offset = round(HYP.latency(i-1))-1; % the offset was the onset of the previous epoch i.e., the offset of the last stage 2 epoch
        thisbout = [thisbout, offset]; %#ok<AGROW> 
        bouts = [bouts; thisbout]; %#ok<AGROW> 
    end
end
% -------------------------------------------------------------------------
% Convert samples to seconds
bouts = bouts./srate;
% -------------------------------------------------------------------------
% Remove segments that are less than 'cutoff' duration
rmidx = diff(bouts') <= cutoff;
bouts(rmidx, :) = [];
end
