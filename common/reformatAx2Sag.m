function [sag] = reformatAx2Sag(ax)
%[sag]=REFORMATAX2SAG(ax) re-arranges axial breast MRI volume into sagittal
%view, without 
%   此处显示详细说明

ax = flip(ax, 3); % head/feet flip

[col,slc,row] = size(ax);
if islogical(ax)
    sag = false(row,col,slc);
else  
    sag = zeros(row,col,slc, 'like',ax);
end

for k = 1:slc
    for j = 1:col
        sag(:,j,k) = ax(j,k,:);
    end
end

end

