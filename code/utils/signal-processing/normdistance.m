function d = normdistance(chanlocs, badchans, i)
% Get the distance between the i'th badchannel and all other channels
d = sqrt(...
    ([chanlocs.X] - chanlocs(badchans(i)).X).^2 + ...
    ([chanlocs.Y] - chanlocs(badchans(i)).Y).^2 + ...
    ([chanlocs.Z] - chanlocs(badchans(i)).Z).^2);
% Interpolated value depends on the squared distance
d = 1./(d.^2);
d(badchans) = 0; % exclude bad channels
d = d ./ sum(d); % normalize to sum = 1
d = ascolumn(d);
end