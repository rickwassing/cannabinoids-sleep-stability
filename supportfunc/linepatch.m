function h = linepatch(Ax, XData, YData, varargin)

if length(size(XData)) > 2
    error('XData must be a vector or 2D matrix')
end

if length(size(YData)) > 2
    error('YData must be a vector or 2D matrix')
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

if size(XData, 2) == 1 && size(XData, 2) ~= size(YData, 2)
    XData = repmat(XData, 1, size(YData, 2));
end

if (size(XData, 1) ~= size(YData, 1)) || (size(XData, 2) ~= size(YData, 2))
    error('X and YData must be equal dimensions')
end

XData = [XData; nan(1, size(XData, 2))];
XData = XData(:);

YData = [YData; nan(1, size(YData, 2))];
YData = YData(:);

h = patch(Ax, 'XData', XData, 'YData', YData, 'LineStyle', '-', varargin{:})

end