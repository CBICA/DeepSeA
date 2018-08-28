function [respMap,smooth_slices] = CWL_likelihood_filteringAx3D(slices,CWL_ROI,frontBounds,half_height,dispOpt)
%[respMap,smooth_slices]=CWL_LIKELIHOOD_FILTERINGAX(slices,CWL_ROI,half_height,dispOpt)

slices = mat2gray(slices);
[row,col,slc] = size(slices);

% smoothing valid/given slices one by one
smooth_slices = zeros(row,col,slc);
for i = 1:slc
%     w         = 5;       % bilateral filter half-width
%     sigma     = [3 0.1]; % bilateral filter standard deviations
%     smooth_slices(:,:,i) = bfilter2(slices(:,:,i),w,sigma);
%     smooth_slices(:,:,i) = anisodiff(slices(:,:,i), 30, 20, 0.1, 1);
    smooth_slices(:,:,i) = wiener2(slices(:,:,i),[7 7]);
end
% smooth_slices = smooth3(slices, 'gaussian');

cubic_len = 2*half_height + 1;
num_voxels = row*col*slc;
% tic
num_samples = 2^22;
if num_samples>num_voxels
    num_samples = num_voxels;
end
p = randperm(num_voxels,num_samples); % have to randomly sample when data is too big
[~,centers]=kmeanspp(smooth_slices(p),3);
% kmeans_of_random_samples = toc
centers=sort(centers);
template = [centers(3)*ones(half_height-1,1); (centers(3)+centers(1))/2; centers(1); (centers(2)+centers(1))/2; centers(2)*ones(half_height-1,1)];
template = template - mean(template(:)); % shift mean to 0
template = repmat(template,1,cubic_len,cubic_len);

respMap = imfilter(smooth_slices, template, 'replicate');
% NOTE: [0,1] scaling should only be done after dark preference
% multiplication; confirmed by experiments!
respMap = respMap.*(1-smooth_slices); % favor dark pixels
% scale using responses only within ROI
amin = min(respMap(CWL_ROI));
amax = max(respMap(CWL_ROI));
respMap = mat2gray(respMap,[amin amax]);
respMap(~CWL_ROI) = 0; % null response outside ROI

% visual debug codes
if exist('dispOpt', 'var') && dispOpt
    implay([mat2gray(slices) ones(row,1,slc) respMap]);
end

end

