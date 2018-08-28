function [respMap,smooth_slices] = CWL_likelihood_filtering(slices,CWL_ROI,leftBounds,half_width,visualize)
%[respMap,smooth_slices]=CWL_LIKELIHOOD_FILTERING(slices,CWL_ROI,leftBounds,half_width,visualize)

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

piece_len = 2*half_width + 1;
% template = [ones(1,half_width) -ones(1,half_width+1)];
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
template = repmat(template,piece_len,1);
% % visual display codes that generate figure for paper
% imwrite(template,'generic_template.png');
template = template - mean(template(:)); % shift mean to 0
% % normalize self response to 1; needed in matching filtering ONLY
% factorsquare = sum(template(:).^2);
% template = template/sqrt(factorsquare);
% rightside = [zeros(1,half_width+1) ones(1,half_width)/half_width];
% rightside = repmat(rightside,piece_len,1);

% slice-wise filtering
respMap = zeros(row,col,slc);
% fun = @(x) rawDist(x,template);
for i = 1:slc
    if leftBounds(i)>0 % skip slices that do not contain breast
        ROI = smooth_slices(:,leftBounds(i)-half_width:end,i);
        rawresp = imfilter(ROI, template, 'replicate');
%         rightresp = imfilter(ROI, rightside, 'replicate');
%         rawresp = nlfilter(ROI,[size(template,1) size(template,2)],fun);
%         rawresp = matchingFiltering(ROI, template, 15);
%         rawresp = rawresp.*rightresp;
        respMap(:,leftBounds(i):end,i) = rawresp(:,1+half_width:end);
    end
end
% NOTE: [0,1] scaling should only be done after dark preference
% multiplication; confirmed by experiments!
respMap = respMap .* (1-smooth_slices); % favor dark pixels
% scale using responses only within ROI
respMap = masked_feature_01norm(respMap, CWL_ROI);
respMap(~CWL_ROI) = 0;

% visual debug codes
if exist('visualize', 'var') && visualize
    implay([mat2gray(slices) ones(row,1,slc) respMap]);
end

end

function t=NCC(V,U)
Vvar=V-mean(V(:)); Uvar=U-mean(U(:));
I=(Vvar.*Uvar)/((sqrt(sum(Vvar(:).^2))*sqrt(sum(Uvar(:).^2)))+eps);
t=sum(I(:))/numel(I);
end

function d = rawDist(V,U)
d = sum(abs(V(:)-U(:)));
end

function responsemap = matchingFiltering(img, K, step)
theta = -30:step:30;
responsemapmx = zeros([size(img) numel(theta)]); % pre-allocation of memory
for i = 1:numel(theta)
    K1 = imrotate(K, theta(i), 'bilinear');
    K1 = K1 - mean(K1(:)); % shift mean to 0
    factorsquare = sum(K1(:).^2);
    K1 = K1/sqrt(factorsquare);
    responsemapmx(:,:,i) = imfilter(img, K1, 'replicate');
end

responsemap = max(responsemapmx, [], 3);
end

