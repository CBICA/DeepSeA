function [gloTemp] = globalAdaptSag3d(CWL_map,paddedVol,respMap,half_height,half_width)
%[gloTemp] = GLOBALADAPTSAG3D(CWL_map,paddedVol,respMap,half_height,half_width)
%   此处显示详细说明

gloTemp = zeros(2*half_height+1,2*half_width+1,2*half_height+1);
sum_w = 0;
linInd = find(CWL_map);
[r,c,s] = ind2sub(size(CWL_map), linInd);
for i = 1:numel(linInd)
    % extract local template and add it to the global template
    gloTemp = gloTemp + respMap(r(i),c(i),s(i)) * paddedVol( r(i)+(0:2*half_height), c(i)+(0:2*half_width), s(i)+(0:2*half_height) );
    sum_w = sum_w + respMap(r(i),c(i),s(i));
end
gloTemp = gloTemp/sum_w;

end

