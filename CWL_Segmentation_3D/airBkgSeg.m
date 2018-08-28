function [mask_init,validSlc,validRow] = airBkgSeg(slices,spacing,szThresh,htThresh)
%[mask_init,validSlc,validRow]=AIR_BKG_SEG(slices,spacing,szThresh,htThresh)
% identifies air/background regions and returns a mask excluding them.
%   INPUTS:
%   - szThresh: area size threshold in determining if a connected component
%     belongs to the initial body mask, in cm2
%   - htThresh: height threshold in determining if connected component
%     belongs to the initial body mask, in terms of ratio w.r.t the no. of
%     image rows

% parse inputs
if ~exist('szThresh', 'var')
    szThresh = 2^2;
end
if ~exist('htThresh', 'var')
    htThresh = 1/3;
end

pxsize = prod(spacing(1:2)/10); % pixel size in cm2

% use double-precision data remapped to [0,1] for processing
slices = mat2gray(slices);

[row,col,slc] = size(slices);

fprintf('### Air/Background Removal Module ###\n');

% pure intensity-based thresholding -- uncomment to use
% thresh = graythresh(slices);
% thresh = air_threshold_calc(slices);
% mask_bright = slices>=thresh;

% 2-d slice-wise filtering
averageMat = slices;
entropyMat = slices;
for i = 1:slc
    averageMat(:,:,i) = imfilter(slices(:,:,i), fspecial('gaussian',[3 3], .5));
    entropyMat(:,:,i) = stdfilt(slices(:,:,i), ones(3));
end

% 3-d filtering
% averageMat = smooth3(slices, 'gaussian');
% entropyMat = stdfilt(slices, ones(3,3,3));

% clustering
% [ctr, pstr] = fcm(averageMat(:), 2);
% [ctr, pstr] = fcm([averageMat(:),entropyMat(:)], 4);
% [~,cls] = max(pstr);
% [cls, ctr] = kmeans([averageMat(:),entropyMat(:)], 4, 'MaxIter', 200); 
% [~,cen] = min(ctr(:,1));
[cls, ctr] = kmeanspp([averageMat(:),entropyMat(:)]', 4); 
[~,cen] = min(ctr(1,:));
mask_bright = false(row,col,slc);
mask_bright(:) = ~(cls==cen);

% 3D morphological processing to keep only the largest CC
mask_init = biggest3dCCinMaskInit(mask_bright);
m = mean(slices(mask_init));
sd = std(slices(mask_init));

% iterative 2D+3D morphing
converged = false;
SE = strel('disk',3);
for iter = 1:15
    mask_ctrl = mask_init;
    
    % Slice-wise 2D morphological processing
    for i = 1:slc
        if any(any(mask_init(:,:,i))) % only process slices with suspect breast
            % break narrow connections before morph operations
            mask_slc = imerode(mask_init(:,:,i),SE);
            CC = bwconncomp(mask_slc,4);
            if CC.NumObjects<1 % in case nothing left after erosion
                mask_init(:,:,i) = false;
                continue;
            end
            
            numPixels = cellfun(@numel, CC.PixelIdxList);
            [numPixels, idx] = sort(numPixels, 'descend');
            mask_slc(:) = false;
            
            % Keep only the largest two components in every slice; used two (2)
            % instead of one due to potentially large breasts and/or skin folding
            [CCPixelRsub,~] = ind2sub([row col],CC.PixelIdxList{idx(1)});
            if pxsize*numPixels(1)>szThresh && range(CCPixelRsub)/row>htThresh
                mask_slc(CC.PixelIdxList{idx(1)}) = true;
            else
                mask_init(:,:,i) = false;
                continue;
            end
            % comment codes below to keep only the largest CC
            if CC.NumObjects>1 && pxsize*numPixels(2)>szThresh
                if mean(CC.PixelIdxList{idx(2)}) > mean(CC.PixelIdxList{idx(1)})
                    % in case the 2nd largest CC is actually body region
                    mask_slc(CC.PixelIdxList{idx(2)}) = true;
                else
                    mask_cand = false(row,col);
                    mask_cand(CC.PixelIdxList{idx(2)}) = true;
                    slice = slices(:,:,i);
                    if ~isAttachedPect(mask_cand,mask_slc) && mean(slice(mask_cand))>(m-sd/2) % the 2nd largest CC should not be dark ghost
                        mask_slc = mask_slc | mask_cand;
                    end
                end
            end
            mask_slc = imfill(imdilate(mask_slc,SE), 'holes');
%             figure(1);
%             imshow(mask_slc);
%             pause;
            mask_init(:,:,i) = mask_slc;
        end
    end
    
    % Then 3D morphological processing again to keep only the largest CC
    mask_init = biggest3dCCinMaskInit(mask_init);
    
    if isequal(mask_ctrl, mask_init)
        disp(['Converged at Iter ' num2str(iter) '.']);
        converged = true;
        break
    end
end

if ~converged
    warning('The iterative morph process did not converge.');
end

% get indices of valid slices that contain breast
validSlc = transpose( find(any(any(mask_init,1),2)) );

% get row range of mask_init in valid slices
validRow = verticalRange_sag(mask_init(:,:,validSlc));

end

function YN = isAttachedPect(mask_cand,mask_slc)
[r1,c1] = find(mask_slc);
[r2,c2] = find(mask_cand);
% Is the candidate shorter than the main component?
if (max(r2)-min(r2))>=(max(r1)-min(r1))
    YN = false;
    return;
end

% Is the candidate strictly to the right of the main component?
c1 = c1(r1<=max(r2)&r1>=min(r1));
if mean(c2)<=mean(c1)
    YN = false;
else
    YN = true;
end

end

