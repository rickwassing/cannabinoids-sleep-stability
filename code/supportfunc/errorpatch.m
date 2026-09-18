function h = errorpatch(Ax, XData, YData, EData, varargin)

if ~isvector(XData)
    error('XData must be a vector')
end

if ~isvector(YData)
    error('YData must be a vector')
end

if ~isvector(EData)
    error('EData must be a vector')
end

if any(size(XData) == 1)
    if size(XData, 1) == 1
        XData = XData';
    end
end

if any(size(YData) == 1)
    if size(YData, 1) == 1
        YData = YData';
    end
end

if any(size(EData) == 1)
    if size(EData, 1) == 1
        EData = EData';
    end
end

if (size(XData, 1) ~= size(YData, 1)) || (size(XData, 2) ~= size(YData, 2))
    error('X and YData must be equal dimensions')
end
XData = [XData; flipud(XData)];
EData = [YData-EData; flipud(YData+EData)];

h = patch(Ax, 'XData', XData, 'YData', EData, 'LineStyle', 'none', varargin{:});

end