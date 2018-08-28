function [cleanmask] = CCFilterRemoveSmallOnes(initmask, N)
% [cleanmask] = CCFilterRemoveSmallOnes(initmask, N)
%   connected component number analysis

CC = bwconncomp(initmask);
numPixels = cellfun(@numel, CC.PixelIdxList);

idx = find(numPixels<=N);
cleanmask = initmask;
for i = 1:length(idx)
    cleanmask(CC.PixelIdxList{idx(i)}) = 0;
end

end

