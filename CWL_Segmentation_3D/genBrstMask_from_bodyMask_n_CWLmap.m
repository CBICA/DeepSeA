function [mask_breast] = genBrstMask_from_bodyMask_n_CWLmap(mask_body,CWL_map)
%[mask_breast]=genBrstMask_from_bodyMask_n_CWLmap(mask_body,CWL_map)
%   此处显示详细说明

mask_breast = mask_body;

for k = 1:size(CWL_map,3)
    [r,c]=find(CWL_map(:,:,k));
    if isempty(r) % null breast mask completely for slices without CWL being detected
        mask_breast(:,:,k) = false;
    else
        if numel(unique(r)) ~= numel(r)
            error('Each row should have only one CWL location!');
        end
        for i = 1:numel(r)
            bodyEnd = find(mask_breast(r(i),:,k),1,'last');
            % Note the dark CWL itself is excluded from the breast mask
            if bodyEnd>c(i) % need to exclude body behind CWL
                mask_breast(r(i),c(i):end,k) = false;
            else % need to include breast in front of CWL
                mask_breast(r(i),bodyEnd+1:c(i)-1,k) = true;
            end
        end
        % slicewise morphological cleaning
        mask_breast(:,:,k) = CCFilterKeepLargestN(mask_breast(:,:,k), 2);
        mask_breast(:,:,k) = CCFilterRemoveSmallOnes(mask_breast(:,:,k), 32^2);
    end
end

% volumetric morphological cleaning
mask_breast = biggestCC(mask_breast);

% some morphological post-processing
mask_breast = imclose(mask_breast, true(3,3,3));
mask_breast = imfill(mask_breast, 'holes');

% for data structure consistency with the orginal by Shandong/Meng-Kang
mask_breast = uint8(mask_breast)*128;

end

