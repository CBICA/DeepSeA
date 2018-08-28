function [mask] = genBrstMask_from_bodyMask_n_CWLmapAx(mask,CWL_map,validCol1,validCol2)
%[mask_breast]=genBrstMask_from_bodyMask_n_CWLmapAx(mask_body,CWL_map,validCol1,validCol2)
%   此处显示详细说明

% clear voxels outside given lateral boundaries of breasts
col2clear = setdiff(1:size(mask,2), [validCol1 validCol2]);
mask(:,col2clear,:) = 0;

for k = 1:size(CWL_map,3)
    [r,c]=find(CWL_map(:,:,k));
    if isempty(c) % null breast mask completely for slices without CWL being detected
        mask(:,:,k) = false;
    else
        if numel(unique(c)) ~= numel(c)
            error('Each column should have only one CWL location!');
        end
        for i = 1:numel(c)
            bodyEnd = find(mask(:,c(i),k),1,'last');
            % Note the dark CWL itself is excluded from the breast mask
            if bodyEnd>r(i) % need to exclude body behind CWL
                mask(r(i):end,c(i),k) = false;
            else % need to include breast in front of CWL
                mask(bodyEnd+1:r(i)-1,c(i),k) = true;
            end
        end
    end
end

% some morphological post-processing
mask = imclose(mask, true(3,3,3));
mask = imfill(mask, 'holes');

% for data structure consistency with the orginal by Shandong/Meng-Kang
mask = uint8(mask)*128;

end

