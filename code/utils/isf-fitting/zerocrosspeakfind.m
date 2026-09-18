function [peak_amps, peak_idx] = zerocrosspeakfind(data)

% Find the zero-crossings
findzc = @(d) unique([1, find((d(1:end-1) <= 0 & d(2:end) > 0) | (d(1:end-1) >= 0 & d(2:end) < 0)), length(d)]);
idx_zc = arrayfun(@(ridx) findzc(data(ridx, :)), 1:size(data, 1), 'UniformOutput', false)';

% Init
peak_idx = [];
peak_amps = [];

for ch = 1:size(data, 1)

    % Init
    ch_peak_idx = [];
    ch_peak_amps = [];

    % Find peaks between consecutive zero crossings
    for i = 1:length(idx_zc{ch})-1

        start_idx = idx_zc{ch}(i);
        end_idx = idx_zc{ch}(i+1);

        % Extract segment between zero crossings
        segment = data(ch, :);
        segment = smoothdata(segment, 'movmean', 25); % corresponding to 0.2 Hz
        segment = segment(start_idx:end_idx);

        % Find max in this segment (only consider positive segments)
        [~, idx_absmax] = max(abs(segment));
        if segment(idx_absmax) > 0 && idx_absmax ~= 1 && idx_absmax ~= length(segment)
            [max_val, max_idx_rel] = max(segment);
            max_idx_abs = start_idx + max_idx_rel - 1;
            max_idx_abs(max_idx_abs == length(data(ch, :))) = [];
            % Store peak time and amplitude
            ch_peak_idx = [ch_peak_idx, max_idx_abs];
            ch_peak_amps = [ch_peak_amps, max_val];
        end
    end

    peak_idx{ch, 1} = ch_peak_idx;
    peak_amps{ch, 1} = ch_peak_amps;
end
