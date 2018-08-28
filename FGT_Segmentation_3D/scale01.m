function [featureVector] = scale01(featureVector, low, high)
%[featureVector]=SCALE01(featureVector, low, high)
%   此处显示详细说明

if low == high
    error('Argument LOW must be smaller than argument HIGH!');
end

if min(featureVector(:))<low
    error('Minimum of given features is smaller than predefined lower boundary and will result in negative values.');
end

featureVector = (double(featureVector) - double(low)) / (double(high) - double(low));

end

