function EEG = executeappending(EEG, varargin)

if nargin < 2
    type = 'add';
else
    type = varargin{:};
end

switch type 
    case 'add'
        EEG.data = [EEG.prepend, EEG.data, EEG.append];
        EEG.pnts = size(EEG.data, 2);
        EEG.xmin = 0;
        EEG.xmax = EEG.pnts/EEG.srate;
        EEG.times = linspace(EEG.xmin, EEG.xmax, EEG.pnts);
        for e = 1:length(EEG.event)
            EEG.event(e).latency = EEG.event(e).latency+size(EEG.prepend, 2);
        end
    case 'remove'
        EEG.data = EEG.data(:, size(EEG.prepend, 2)+1:end-size(EEG.append, 2));
        EEG.pnts = size(EEG.data, 2);
        EEG.xmin = 0;
        EEG.xmax = EEG.pnts/EEG.srate;
        EEG.times = linspace(EEG.xmin, EEG.xmax, EEG.pnts);
        for e = 1:length(EEG.event)
            EEG.event(e).latency = EEG.event(e).latency-size(EEG.prepend, 2);
        end
end

end