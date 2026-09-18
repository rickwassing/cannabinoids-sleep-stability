function [sig, cmpx] = predict_isf(sig, const, doPlot)

if nargin < 3
    doPlot = false;
end

N = abs(const.crop_aro(2))-30*sig.srate;

data_pred = nan(sig.nbchan, N);
cmpx_pred = nan(sig.nbchan, N);
cmpx_sig = [];
for chan = 1:sig.nbchan
    [data_pred(chan, :), cmpx_pred(chan, :), cmpx_sig(chan, :)] = ...
        predict_infraslow_eeg(detrend(double(sig.data(chan,:)), 0), ...
        'sampling_rate', sig.srate, ...
        'prediction_length', N, ...
        'do_plot', doPlot, ...
        'padding_samples', length(sig.prepend));
    if doPlot
        drawnow();
    end
end

cmpx = [cmpx_sig, cmpx_pred];

sig = rmappend(sig);

sig.data = [sig.data, data_pred];
sig.pnts = size(sig.data, 2);
sig.xmin = 0;
sig.xmax = (sig.pnts-1)/sig.srate;
sig.times = linspace(sig.xmin, sig.xmax, sig.pnts);

end