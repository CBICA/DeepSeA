function [BW] = biggestCC(BW)
%[BW] = BIGGESTCC(BW)
%   Keep only the biggest connected component in a binary image (2D or 3D)

CC = bwconncomp(BW);
numPixels = cellfun(@numel,CC.PixelIdxList);
[~,idx] = max(numPixels);
BW(:) = false;
BW(CC.PixelIdxList{idx}) = 1;

end

