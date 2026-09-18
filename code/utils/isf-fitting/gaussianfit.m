function [res, offset] = gaussianfit(XData, YData, N)
% -------------------------------------------------------------------------
% Prepare input data, columnize, remove NaN's etc.
XData = ascolumn(XData);
YData = ascolumn(YData);
[xData, yData] = prepareCurveData(XData, YData);
yData(xData < 0.0075) = linspace(median(yData), median(yData), sum(xData < 0.0075));
offset = mean(yData(xData >= 0.059999));
yData = yData - offset;
% -------------------------------------------------------------------------
% Find initial fitting values
[~, peakFreqIdx] = max(yData(xData > 0.0075 & xData <= 0.04));
peakFreqIdx = find(xData > 0.0075, 1, 'first') + peakFreqIdx;
peakFreq = xData(peakFreqIdx);
% Minimum standard deviation is 2 times the frequency resolution
minsd = 2*mean(diff(xData));
% -------------------------------------------------------------------------
% Set up fittype and options.
opts = fitoptions('Method', 'NonlinearLeastSquares');
opts.Display = 'Off';
% -------------------------------------------------------------------------
% Starting points and lower and upper limits of the parameters
% Amplitude, Mu, and SD of Gaus1, and then Amp, Mu and SD of Gaus2, etc.
res = struct();
adjr = [];
for i = 1:N
    switch i
        case 1
            ft = fittype('gauss1');
            opts.StartPoint = [...
                max(yData), peakFreq, 0.02]; % amp, mu, std
            opts.Lower = [...
                0.01*max(yData), 0.0075+minsd, minsd];
            opts.Upper = [...
                2*max(yData), 0.1-minsd, 0.2];
        case 2
            ft = fittype('gauss2');
            opts.StartPoint = [...
                max(yData), peakFreq, 0.02, ...
                0.5*max(yData), 0.06, 0.01];
            opts.Lower = [...
                0.01*max(yData), 0.0075+minsd, minsd, ...
                0.01*max(yData), 0.0267+minsd, minsd];
            opts.Upper = [...
                2*max(yData), 0.0333-minsd, 0.2, ...
                2*max(yData), 0.1000, 0.2];
    end
    % -------------------------------------------------------------------------
    % Fit model to data.
    [res(i).fobj, gof] = fit(xData, yData, ft, opts);
    adjr = [adjr, gof.adjrsquare]; %#ok<AGROW>
end
% -------------------------------------------------------------------------
% Select most parsimonious model
[~, idx] = max(adjr);
res = res(idx).fobj;
