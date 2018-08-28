function [CWLs,CWL_map] = CWL_localiseAx(validSlcID,index1,index2,slices,spacing,CWL_ROI,frontBounds,filter_size_half,lambdaAdjuster,iso_spacing,dimension,visualize)
%CWL_LOCALISEAX Summary of this function goes here
%   Detailed explanation goes here

% parse inputs
if ~exist('visualize', 'var')
    visualize = false;
end

% find the lateral bounds to cut matrices to pass into subsequent function
lateral_lb = min(index1);
lateral_ub = max(index2);

d = dMapAx(slices(:,lateral_lb:lateral_ub,validSlcID),spacing, ...
    CWL_ROI(:,lateral_lb:lateral_ub,validSlcID),frontBounds(validSlcID), ...
    filter_size_half,lambdaAdjuster,iso_spacing,dimension,visualize);

CWLs = cell(size(slices,3),1);
for i = 1:size(d,2)
    CWLs{validSlcID(i)} = [transpose(lateral_lb:lateral_ub) d(:,i)];
end

[row,col,slc] = size(slices);
CWL_map = false(row,col,slc);
[slcSub,colSub] = meshgrid(validSlcID,lateral_lb:lateral_ub);
CWL_map(sub2ind([row,col,slc],round(d(:)),colSub(:),slcSub(:)))...
    = true;

% visual debug codes
if exist('visualize','var') && visualize
    slices_4disp = mat2gray(slices);
    slices_4disp(CWL_map) = 1;
    implay(slices_4disp);
end

end

