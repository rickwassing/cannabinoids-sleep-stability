function event = storeoriglatency(event)
for i = 1:length(event)
    event(i).latency = round(event(i).latency);
    event(i).origlatency = event(i).latency;
end
end
