function [VT,FGT,VBD] = calMetrics(mask1,spacing)
%[VT,FGT,VBD] = CALMETRICS(mask1,spacing)
%   VT:  whole breast volume (cm^3)
%   FGT: absolute volume (cm^3) of breast dense tissue
%   VBD: percent breast density (%)

voxel = prod(spacing/10);

if any(mask1(:))
    VT = sum(mask1(:)>0) * voxel;
else % mask1 has not been initialilzed at all
    VT = 0;
    FGT = 0;
    VBD = 0;
    return
end

if any(mask1(:)==255)
    FGT = sum(mask1(:)==255) * voxel;
    VBD = FGT/VT *100;
else
    FGT = 0;
    VBD = 0;
end

end

