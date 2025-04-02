function EEG = convspindles(EEG, s)

x = -6*s:1/EEG.srate:6*s;
g = exp(-1.*((x.^2)./(2*s.^2)));

d_ = nan(EEG.nbchan, EEG.pnts+length(g)-1, 'single');
rt = now();
for i = 1:EEG.nbchan
    d_(i, :) = conv(EEG.data(i, :), g);
    rt = remainingTime(rt, EEG.nbchan);
end

EEG.data = d_(:, (length(g)-1)/2:end-(length(g)-1)/2-1);

end