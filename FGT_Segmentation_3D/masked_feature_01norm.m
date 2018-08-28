function [featureMap] = masked_feature_01norm(featureMap, mask)
%[featureMap] = masked_feature_01norm(featureMap, mask)
%   此处显示详细说明

low = min( featureMap(mask) );
high = max( featureMap(mask) );
featureMap(mask) = scale01(featureMap(mask), low, high);

end

