function sig = rmappend(sig)

npre = size(sig.prepend, 2);
npost = size(sig.prepend, 2);

sig.data = sig.data(:, npre+1:end-npost);
sig.pnts = size(sig.data, 2);
sig.xmin = 0;
sig.xmax = (sig.pnts-1)/sig.srate;
sig.times = linspace(sig.xmin, sig.xmax, sig.pnts);
sig = rmfield(sig, 'append');
sig = rmfield(sig, 'prepend');
end