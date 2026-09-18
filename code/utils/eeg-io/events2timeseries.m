function x = events2timeseries(event, xmin, xmax, srate)

times = xmin:1/srate:xmax;
x = false(1, length(times));
for i = 1:length(event)
    x(round(event(i).latency):round(event(i).latency + event(i).duration)) = true;
end
x = x(1:length(times));
end
