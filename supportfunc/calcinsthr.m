function EEG = calcinsthr(EEG)
% -------------------------------------------------------------------------
% Find ECG channel
idx_ecg = find(strcmpi({EEG.chanlocs.type}, 'ecg'));
% -------------------------------------------------------------------------
% Find QRS complexes
[~, idx_qrs, ~, ecg, Fig] = ecg2hr_pantompkin(EEG.data(idx_ecg, :), EEG.srate, 1); % 1 = do plot
savefig(Fig, sprintf('inspect/ecg/%s.fig', EEG.setname))
% -------------------------------------------------------------------------
% Convert to BPM
bpm = 60 ./ (diff(idx_qrs) ./ EEG.srate);
% -------------------------------------------------------------------------
% Cleaning
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Remove -20 and +200 BPM
bpm(bpm < 20) = NaN;
bpm(bpm > 200) = NaN;
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Calculate beat-to-beat changes in HR (percentage)
for iteration = 1:3
    db2b = 0;
    for i = 2:length(bpm)
        db2b(i) = (bpm(i) - bpm(i-1)) / bpm(i-1);
    end
    % Z-score, and find outliers (3 standard deviations)
    idx_rm = abs(nanzscore(db2b)) > 3;
    % Remove beat-to-beat outliers
    bpm(idx_rm) = NaN;
    % Interpolate missing values
    idx_nan = find(isnan(bpm));
    if isempty(idx_nan)
        break
    end
    idx_nan = [idx_nan(1), idx_nan(find(diff(idx_nan) ~= 1)+1)];
    for i = 1:length(idx_nan)
        % Find the next non-nan value
        n = find(~isnan(bpm(idx_nan(i):end)), 1, 'first') - 1;
        if isempty(n) % must be end of the signal
            bpm(idx_nan(i):end) = median(bpm, 'omitnan');
            continue
        end
        % Extract the last-known bpm
        if idx_nan(i) == 1
            x1 = median(bpm, 'omitnan');
        else
            x1 = bpm(idx_nan(i)-1);
        end
        % Extract the first-known bpm following the missing data
        if idx_nan(i) == length(bpm)
            x2 = median(bpm, 'omitnan');
        else
            x2 = bpm(idx_nan(i)+n);
        end
        x = linspace(x1, x2, n+2);
        x([1, end]) = [];
        bpm(idx_nan(i):idx_nan(i)+n-1) = x;
    end
end
% -------------------------------------------------------------------------
% Construct instantaneous HR timeseries
hr = nan(1, EEG.pnts, 'single');
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Assume the HR for the first detected QRS
hr(1:idx_qrs(1)) = bpm(1);
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% For each detected QRS, create the instantaneous HR
for i = 1:length(idx_qrs)
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Start index
    x1 = idx_qrs(i);
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % End index
    if (i+1) > length(idx_qrs)
        x2 = EEG.pnts;
    else
        x2 = idx_qrs(i+1);
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Start value
    if i == 1
        y1 = bpm(1);
    else
        y1 = bpm(i-1);
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % End value
    if i > length(bpm)
        y2 = bpm(end);
    else
        y2 = bpm(i);
    end
    % - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    % Insert into timeseries vector
    if (x2 - x1 + 1) > 3*EEG.srate
        hr(x1:x2) = nan;
    else    
        hr(x1:x2) = linspace(y1, y2, x2-x1+1);
    end
end
% -------------------------------------------------------------------------
% Append to data
EEG.data = [EEG.data; hr];
EEG.nbchan = EEG.nbchan+1;
EEG.chanlocs(end+1).labels = 'HR';
EEG.chanlocs(end).ref = 'average';
EEG.chanlocs(end).type = 'HR';
EEG.chanlocs(end).unit = 'bpm';

end