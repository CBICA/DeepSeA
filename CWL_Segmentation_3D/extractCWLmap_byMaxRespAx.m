function [true_positive] = extractCWLmap_byMaxRespAx(respMap,slices)
%[CWL_map]=EXTRACTCWLMAP_BYMAXRESP(respMap,slices)
%   此处显示详细说明

[row,col,slc] = size(respMap);

[~,d0] = max(respMap,[],1);
d0 = reshape(d0,col,slc);
if any(d0(:)==1)
    warning('Found all equal row(s) in response map!')
%     [c,s,v] = find(d0~=1);
%     [Cq,Sq] = meshgrid(1:col,1:slc);
%     F = scatteredInterpolant(c,s,v,'nearest','nearest');
%     d0 = F(Cq,Sq);
    d0(d0==1) = row;
end

CWL_map = false(row,col,slc);
[slcSub,colSub] = meshgrid(1:slc,1:col);
CWL_map(sub2ind([row,col,slc],d0(:),colSub(:),slcSub(:)))...
    = true;

% false positive removal: keep only the largest CC
true_positive = biggestCC(CWL_map);

% visual debug codes
if exist('slices', 'var') % used as a flag for visual debugging
    false_positive = CWL_map & ~true_positive;
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

