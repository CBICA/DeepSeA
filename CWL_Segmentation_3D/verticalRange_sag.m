function [validRow] = verticalRange_sag(mask3D)
%[validRow] = VERTICALRANGE_SAG(mask3D)
%   此处显示详细说明

linearInd = find(mask3D);
[rowSub,~,~] = ind2sub(size(mask3D), linearInd);
validRow = min(rowSub):max(rowSub);

end

