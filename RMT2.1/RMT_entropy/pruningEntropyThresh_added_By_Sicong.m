function [features,dpscale,To_Store] = pruningEntropyThresh_added_By_Sicong(features, dpscale, percentile, data)
% prune feature based on the percentile of max entropy value
time_scope = features(4, :) * 3; 
Feature_Entropy = zeros(1, size(features, 2));
for i = 1 : size(features, 2)
    Candidate_Entropic_Feature = features(:,i);
    Candidate_DepdScale = dpscale(:,i);
    Candidate_Time_Range = (round((Candidate_Entropic_Feature(2, 1) - time_scope(1, i))) : (round((Candidate_Entropic_Feature(2, 1) + time_scope(1, i)))));
    Quantized_Data = globalQuantization(data(Candidate_DepdScale(Candidate_DepdScale > 0, 1), Candidate_Time_Range((Candidate_Time_Range > 0 & Candidate_Time_Range <= size(data, 2)))));
    EntropyF1 = EntropySingVariate_mex(Quantized_Data', -Inf);
    Feature_Entropy(i) = EntropyF1;
end
maximum = max(Feature_Entropy);
thr = maximum * percentile;
To_Remove = Feature_Entropy >= thr;
To_Store = ones(1, size(features, 2)) - To_Remove;
features = features(:, Feature_Entropy <= thr);
dpscale = dpscale(:, Feature_Entropy <= thr);
end