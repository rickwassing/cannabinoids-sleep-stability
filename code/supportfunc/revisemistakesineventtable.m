function t = revisemistakesineventtable(kv, t)
if strcmpi(kv.sub, 'r0072ed')
    t(end, :) = [];
    t.trial_type{24} = 'loff';
    t.trial_type{1215} = 'lon';
end
if strcmpi(kv.sub, 'r0081sr')
    t.trial_type{32} = 'loff';
    t.trial_type{1300} = 'lon';
end
if strcmpi(kv.sub, 'r0082sr')
    t.trial_type{1027} = 'wake';
end
if strcmpi(kv.sub, 'r0092dm')
    t.trial_type{29} = 'loff';
    t.trial_type{1359} = 'lon';
end
if strcmpi(kv.sub, 'r0101mc')
    t.trial_type{27} = 'loff';
end
if strcmpi(kv.sub, 'r0102mc')
    t.trial_type{25} = 'loff';
end
if strcmpi(kv.sub, 'r0111vn')
    t.trial_type{25} = 'loff';
    t.trial_type{1455} = 'lon';
end
if strcmpi(kv.sub, 'r0121lg')
    t.trial_type{34} = 'loff';
end
if strcmpi(kv.sub, 'r0132mg')
    t.trial_type{39} = 'loff';
end
if strcmpi(kv.sub, 'r0141ai')
    t.trial_type{27} = 'loff';
    t.trial_type{1426} = 'lon';
end
if strcmpi(kv.sub, 'r0142ai')
    t.trial_type{112} = 'wake';
    t.trial_type{1446} = 'lon';
end
if strcmpi(kv.sub, 'r0151sp')
    t.trial_type{31} = 'loff';
    t.trial_type{1393} = 'lon';
end
if strcmpi(kv.sub, 'r0161ef')
    t.trial_type{24} = 'loff';
end
if strcmpi(kv.sub, 'r0162ef')
    t.trial_type{31} = 'loff';
end
if strcmpi(kv.sub, 'r0171pm')
    t.trial_type{24} = 'loff';
end
if strcmpi(kv.sub, 'r0172pm')
    t.trial_type{1753} = 'wake';
    t.trial_type{1759} = 'wake';
end
if strcmpi(kv.sub, 'r0182dr')
    t(1190, :) = [];
end
if strcmpi(kv.sub, 'r0191kj')
    t.trial_type{29} = 'loff';
end
if strcmpi(kv.sub, 'r0201mj')
    t.trial_type{1353} = 'lon';
end
% Check that we only have one 'loff' and one 'lon'
if sum(strcmpi(t.trial_type, 'loff')) > 1 || sum(strcmpi(t.trial_type, 'lon')) > 1
    error('More than one lights off/on event')
end

end