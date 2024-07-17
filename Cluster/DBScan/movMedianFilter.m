function [medians] = movMedianFilter(data,windowSize)
paddata = padarray(data,windowSize-1,'both','replicate');
medians= [];

for i=1:length(paddata)-windowSize
    medians=[medians,median(paddata(i:i+windowSize))];
end
medians=medians(windowSize-1:end);