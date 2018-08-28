function [CWLs] = cutCWL2brst(CWLs,mask_body,validSlc)
%[CWLs]=CUTCWL2BRST(CWLs,mask_body,validSlc)
%   此处显示详细说明

for k = validSlc
    [rowSub,~] = find(mask_body(:,:,k));
    top = min(rowSub);
    bttm = max(rowSub);
    
    CWLs{k}(CWLs{k}(:,2)<top,:) = [];
    CWLs{k}(CWLs{k}(:,2)>bttm,:) = [];
end

end

