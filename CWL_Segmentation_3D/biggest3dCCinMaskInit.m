function [BW] = biggest3dCCinMaskInit(BW)
%BIGGEST3DCCINMASKINIT 3D morphological processing which keeps only the 
%biggest CC in BW
%   Detailed explanation goes here

SE = true(3,3,3);

BW = imerode(BW,SE); % break narrow connections before morph operations
CC = bwconncomp(BW); % find connected components

% Find and keep only the largest component
numPixels = cellfun(@numel, CC.PixelIdxList);
[~, idx] = max(numPixels);
BW(:) = false;
BW(CC.PixelIdxList{idx}) = true;

BW = imdilate(BW,SE); % restore the 3D volume

end

