function q = doubleflip(p, type)

if ~any(size(p) == 1)
    if size(p, 1) < size(p, 2)
        warning('It seems your timeseries data is not one trace per column!')
    end
end

if size(p, 1) == 1
    p = p';
end

q = -1.*flip(p, 1);

switch type
    case 'prepend'
        q = q + repmat(2*p(1, :), size(p, 1), 1);
    case 'append'
        q = q + repmat(2*p(end, :), size(p, 1), 1);
    case 'none'
        % do nothing
end

end