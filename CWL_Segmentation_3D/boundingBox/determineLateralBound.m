function [validSlc,validRow] = determineLateralBound(laterality,portion,masks,slices,spacing,validSlc,validRow,visualize)
%[validSlc,validRow]=DETERMINEOUTERBOUND(laterality,portion,masks,slices,spacing,validSlc,validRow[,visualize=false])
%   Detailed explanation goes here

if ~exist('visualize','var')
    visualize = false;
end

% call the mutual core function shared from axial view, using arguments
% re-arranged from sagittal to axial view (slices, masks and valid rows)
validRow = size(masks, 1) - validRow(end:-1:1) + 1;
[validSlc,validRow] = determineLateralBoundCore(laterality,portion,reformatSag2Ax(masks>0),...
    reformatSag2Ax(slices),[spacing(2),spacing(3),spacing(1)],validSlc,validRow,visualize);
validRow = size(masks, 1) - validRow(end:-1:1) + 1;

end

