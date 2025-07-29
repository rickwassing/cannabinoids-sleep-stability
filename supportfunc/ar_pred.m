function x_pred = ar_pred(x, k, order)

if nargin < 3
    order = k;
end

% Recursive function to process matrix data (channels must be columns)
if ~any(size(x) == 1)
    if size(x, 1) < size(x, 2)
        warning('It seems your signals are shorter than the number of channels. Is that right?')
    end
    x_pred = nan(k, size(x, 2));
    for ch = 1:size(x, 2)
        x_pred(:, ch) = ar_pred(x(:, ch), k, order);
    end
    return
end

if size(x, 1) == 1
    x = x';
end
x = double(x);

N = length(x);

% Estimate AR coefficients using Yule-Walker
[r, lags] = xcorr(x, order, 'biased');  % Autocorrelation
R = toeplitz(r(order+1:end-1));         % Autocorrelation matrix
rhs = r(order+2:end);                   % Right-hand side
a = R \ rhs;                        % AR coefficients (without leading 1)
a = [1; -a];                        % Include leading 1 for filter convention

% Forecast forward by k steps
x_pred = zeros(k,1);
x_ext = [x; x_pred];  % Temporary placeholder

for i = 1:k
    x_ext(N+i) = -a(2:end)' * x_ext(N+i-1:-1:N+i-order);  % Predict using past p samples
end

x_pred = x_ext(end-k+1:end);

end