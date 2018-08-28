function [respMap,smooth_slices] = CWL_likelihood_filtering3D(slices,CWL_ROI,leftBounds,half_width,visualize)
%[respMap,smooth_slices]=CWL_LIKELIHOOD_FILTERING3D(slices,CWL_ROI,leftBounds,half_width,visualize)

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

cubic_len = 2*half_width + 1;
n = row*col*slc;
% tic
k = 2^22;
if k>n
    k = n;
end
p = randperm(n,k); % have to do random sampling when data is too big
[~,C]=kmeanspp(smooth_slices(p),3);
% kmeans_of_random_samples = toc
C=sort(C);
template = [C(3)*ones(1,half_width-1), (C(3)+C(1))/2, C(1), (C(2)+C(1))/2, C(2)*ones(1,half_width-1)];
template = template - mean(template(:)); % shift mean to 0
template = repmat(template,cubic_len,1,cubic_len);

respMap = imfilter(smooth_slices, template, 'replicate');
% NOTE: [0,1] scaling should only be done after dark preference
% multiplication; confirmed by experiments!
respMap = respMap.*(1-smooth_slices); % favor dark pixels
% scale using responses only within ROI
respMap = masked_feature_01norm(respMap, CWL_ROI);
respMap(~CWL_ROI) = 0; % null response outside ROI

% visual debug codes
if exist('visualize', 'var') && visualize
    implay([mat2gray(slices) ones(row,1,slc) respMap]);
end

end

