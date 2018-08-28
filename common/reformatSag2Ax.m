function [ax] = reformatSag2Ax(sag)
%[ax]=REFORMATSAG2AX(sag)
%   Detailed explanation goes here

sag = flipud(sag); % head/feet flip

[slc,row,col] = size(sag);
if islogical(sag)
    ax = false(row,col,slc);
else
    ax = zeros(row,col,slc, 'like', sag);
end

for k = 1:slc
    for j = 1:col
        ax(:,j,k) = sag(k,:,j);
    end
end

end

