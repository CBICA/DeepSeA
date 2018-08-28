function [cleanmask] = CCFilterKeepLargestN(initmask, N)
% [cleanmask] = CCFilterKeepLargestN(initmask, N)
%   Keeps only largest N connected components

CC = bwconncomp(initmask);
numPixels = cellfun(@numel, CC.PixelIdxList);

% keep only the N largest components
[~, idx] = sort(numPixels, 'descend');
cleanmask = false(size(initmask));
if numel(idx) <= N
    cleanmask = initmask;
else
    for i = 1:N
        cleanmask(CC.PixelIdxList{idx(i)}) = 1;
    end
end

end

