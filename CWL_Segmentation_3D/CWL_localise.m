function [CWLs,CWL_map] = CWL_localise(validSlcID,validRowID,slices,spacing,CWL_ROI,leftBounds,filter_size_half,lambdaAdjuster,iso_spacing,Tdimension,visualize)
%[CWLs,CWL_map]=CWL_LOCALISE(validSlcID,validRowID,slices,spacing,CWL_ROI,leftBounds,filter_size_half,lambdaAdjuster,iso_spacing,visualize)
%   此处显示详细说明

% parse inputs
if ~exist('visualize', 'var')
    visualize = false;
end

[row,col,slc] = size(slices);

% the main functioning function
d = dMap(slices(validRowID,:,validSlcID),spacing,CWL_ROI(validRowID,:,validSlcID),...
    leftBounds(validSlcID),filter_size_half,lambdaAdjuster,iso_spacing,Tdimension,visualize);

% extract CWLs from depth map
CWLs = cell(size(slices,3),1);
for i = 1:size(d,2)
    CWLs{validSlcID(i)} = [d(:,i) transpose(validRowID)];
end

% construct CWL map from depth map
CWL_map = false(row,col,slc);
[slcSub,rowSub] = meshgrid(validSlcID,validRowID);
CWL_map(sub2ind([row,col,slc],rowSub(:),round(d(:)),slcSub(:)))...
    = true;

% visual debug codes
if exist('visualize', 'var') && visualize
    slices_4disp = mat2gray(slices);
    slices_4disp(CWL_map) = 1;
    implay(slices_4disp);
end

end

