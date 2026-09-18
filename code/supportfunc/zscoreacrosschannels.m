function data = zscoreacrosschannels(data)
if size(data, 2) < size(data, 1)
    data = data'; % make sure the columns are channels
end
mu = data(:);
sd = std(mu);
mu = mean(mu);
data = (data-mu)./sd;

end