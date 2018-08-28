function [true_positive] = extractCWLmap_byMaxResp(respMap,slices)
%[CWL_map]=EXTRACTCWLMAP_BYMAXRESP(respMap,slices)
%   此处显示详细说明

[row,col,slc] = size(respMap);

[~,d0] = max(respMap,[],2);
if any(d0(:)==1)
    warning('Found all-1 row(s) in response map!')
    d0(d0==1) = col;
end
d0 = reshape(d0,row,slc);

CWL_map = false(row,col,slc);
[slcSub,rowSub] = meshgrid(1:slc,1:row);
CWL_map(sub2ind([row,col,slc],rowSub(:),round(d0(:)),slcSub(:)))...
    = true;

% false d0 removal: keep only the largest CC
true_positive = biggestCC(CWL_map);

% visual debug codes
if exist('slices', 'var') % used as a flag for visual debugging
    false_positive = CWL_map & ~true_positive;
%     % Uncomment codes below to make the display thicker, if needed
%     true_positive_bold = imdilate(true_positive,[1 1 0]);
%     false_positive_bold = imdilate(false_positive,[1 1 0]);
    overlaid = zeros(row,col,3,slc,'uint8');
    for i = 1:slc
        in_uint8 = im2uint8(mat2gray(slices(:,:,i)));
        tmpslc = imoverlayMathWorks(in_uint8, ...
            true_positive(:,:,i), [0 1 0]);
        %     tmpslc = imoverlayMathWorks(tmpslc, ...
        %         bwperim(bwmorph(mask_body(:,:,i),'open')), [1 1 0]);
        overlaid(:,:,:,i) = imoverlayMathWorks(tmpslc, ...
            false_positive(:,:,i), [1 0 0]);
    end
    implay(overlaid);
end

end

