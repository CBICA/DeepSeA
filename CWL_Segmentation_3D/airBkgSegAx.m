function [mask_init,validSlc,index1,index2] = airBkgSegAx(slices,rightInd,leftInd,spacing,AThresh,N)
%[mask_init,validSlc,index1,index2] = air_bkg_segAx(slices,rightInd,leftInd,spacing,AThresh,N)

pixelA = prod(spacing(1:2)/10); % pixel size in cm2

% parse inputs
if ~exist('AThresh','var')
    AThresh = 2.5^2; % default area threashold in cm2
end
if ~exist('N','var')
    N = 3;
end

% use double-precision data remapped to [0,1] for processing
slices = normalize01(slices);

[row,col,slc] = size(slices);

fprintf('### Air/Background Removal Module ###\n');

% pure intensity-based thresholding -- uncomment to use
% thresh = graythresh(slices);
% thresh = air_threshold_calc(slices);
% mask_bright = slices>=thresh;

% clustering
% % old codes below for slice-wise filtering
% smoothMat = slices;
% stdMat = slices;
% for i = 1:slc
%     smoothMat(:,:,i) = imfilter(slices(:,:,i), fspecial('gaussian',[3 3], .5));
%     stdMat(:,:,i) = stdfilt(slices(:,:,i), ones(3));
% end

% [ctr, pstr] = fcm(averageMat(:), 2);
% [ctr, pstr] = fcm([averageMat(:),entropyMat(:)], 4);
% [~,cls] = max(pstr);
% [cls, ctr] = kmeans([averageMat(:),entropyMat(:)], 4, 'MaxIter', 200); 
% [~,cen] = min(ctr(:,1));

% tried to add [x y z] coordinates for clustering; however, it didn't work
% out
% [X,Y,Z] = meshgrid(1:col, 1:row, 1:slc);
% X = X/max(X(:)) * 1e-1;
% Y = Y/max(Y(:)) * 1e-1;
% Z = Z/max(Z(:)) * 1e-1;

smoothMat = smooth3(slices, 'gaussian');
stdMat = stdfilt(slices, ones(3,3,3));
% [cls, ctr] = kmeanspp(smoothMat(:)', 2); % pure intensity based clustering
[cls, ctr] = kmeanspp([smoothMat(:),stdMat(:)]', 4); 
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
for iter = 1:12
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
            
            % Keep only the largest three components in every slice; used 3
            % instead of one due to potentially large breasts and/or skin folding
%             [CCPixelRsub,~] = ind2sub([row col],CC.PixelIdxList{idx(1)});
%             htThresh = .5;
            if pixelA*numPixels(1)>AThresh %&& range(CCPixelRsub)/row>htThresh
                mask_slc(CC.PixelIdxList{idx(1)}) = true;
            else
                mask_init(:,:,i) = false;
                continue;
            end
            N_i = min(CC.NumObjects, N);
            for n = 2:N_i
                mask_major = mask_slc;
                if pixelA*numPixels(n)>AThresh
                    mask_cand = false(row,col);
                    mask_cand(CC.PixelIdxList{idx(n)}) = true;
                    slice = slices(:,:,i);
                    if isBrst(mask_cand, mask_major) && mean(slice(mask_cand))>(m-2*sd) % the 2nd & 3rd largest CC should not be dark ghost
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

% find the lateral bounds containing breasts
linearInd = find(mask_init(:,:,validSlc));
[~,subCol,~] = ind2sub(size(mask_init(:,:,validSlc)), linearInd);
subCol = transpose(unique(subCol));
index1 = subCol(subCol<=max(rightInd));
index2 = subCol(subCol>=min(leftInd));

end

function [true_or_false]=isBrst(mask_cand, mask_main)
% a simple function that judges if a small CC belongs to breast based on
% relative position
[r_cand,c_cand] = find(mask_cand);
[r_main,c_main] = find(mask_main);

[r_cand_min, i_cand_min] = min(r_cand);
c_cand_min = c_cand(i_cand_min);
r_main_col = r_main(c_main==c_cand_min);
if ~isempty(r_main_col) && r_cand_min > max(r_main_col)
    true_or_false = false;
else
    true_or_false = true;
end

end

