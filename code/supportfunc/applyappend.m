function sig = applyappend(sig)

sig.data = [sig.prepend, sig.data, sig.append];
sig.pnts = size(sig.data, 2);
sig.xmin = 0;
sig.xmax = (sig.pnts-1)/sig.srate;
sig.times = linspace(sig.xmin, sig.xmax, sig.pnts);

end