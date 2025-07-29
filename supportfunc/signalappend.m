function sig = signalappend(sig, append_type, ford_mult, filtcfg)
switch append_type
    case 'flip'
        % Append a flipped copy of itself at the end of the signal
        if sig.pnts < ford_mult*filtcfg.order
            sig.prepend = zeros(sig.nbchan, ford_mult*filtcfg.order);
            sig.append = zeros(sig.nbchan, ford_mult*filtcfg.order);
            sig.prepend(:, (end-sig.pnts)+1:end) = doubleflip(sig.data', 'none')'; %-1.*flip(sig.data, 2);
            sig.append(:, 1:sig.pnts) = doubleflip(sig.data', 'none')'; %-1.*flip(sig.data, 2);
        else
            sig.prepend = doubleflip(sig.data(:, 1:ford_mult*filtcfg.order)', 'none');
            sig.append = doubleflip(sig.data(:, end-ford_mult*filtcfg.order+1:end)', 'none');
        end
    case 'zeros'
        % Append with zeros
        sig.prepend = zeros(sig.nbchan, ford_mult*filtcfg.order);
        sig.append = zeros(sig.nbchan, ford_mult*filtcfg.order);
    case 'autoreg'
        % Autoregressive prediction
        if sig.pnts < ford_mult*filtcfg.order
            sig.prepend = flip(ar_pred(flip(sig.data, 2)', ford_mult*filtcfg.order, sig.pnts)', 2);
            sig.append = ar_pred(sig.data', ford_mult*filtcfg.order, sig.pnts)';
        else
            sig.prepend = flip(ar_pred(flip(sig.data, 2)', ford_mult*filtcfg.order)', 2);
            sig.append = ar_pred(sig.data', ford_mult*filtcfg.order)';
        end
    case 'none'
        sig.prepend = [];
        sig.append = [];
    otherwise
        error('Specify an append type.')
end
end