function [d_field] = CWL_mapTo_d_field(CWL_map)
%[d_field]=CWL_MAPTO_D_FIELD(CWL_map) 此处显示有关此函数的摘要
%   此处显示详细说明

[row,~,slc] = size(CWL_map);

d_field = zeros(row, slc);
for k = 1:slc
    [r,c]=find(CWL_map(:,:,k));
    d_field(r,k)=c;
end

end

