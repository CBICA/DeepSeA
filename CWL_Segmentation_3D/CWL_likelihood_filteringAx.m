function [respMap,smooth_slices] = CWL_likelihood_filteringAx(slices,CWL_ROI,frontBounds,half_height,dispOpt)
%[respMap,smooth_slices]=CWL_LIKELIHOOD_FILTERINGAX(slices,CWL_ROI,leftBounds,half_height,dispOpt)

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

piece_width = 2*half_height + 1;
% template = [ones(1,half_height) -ones(1,half_height+1)];
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
template = [C(3)*ones(half_height-1,1); (C(3)+C(1))/2; C(1); (C(2)+C(1))/2; C(2)*ones(half_height-1,1)];
template = template - mean(template(:)); % shift mean to 0
template = repmat(template,1,piece_width);
% % normalize self response to 1; needed in matching filtering ONLY
% factorsquare = sum(template(:).^2);
% template = template/sqrt(factorsquare);
% rightside = [zeros(1,half_width+1) ones(1,half_width)/half_width];
% rightside = repmat(rightside,piece_len,1);

respMap = zeros(row,col,slc);
% fun = @(x) rawDist(x,template);
for i = 1:slc
    if frontBounds(i)>0 % only filter slices that contain breast
        ROI = smooth_slices(frontBounds(i)-half_height:end,:,i);
        rawresp = imfilter(ROI, template, 'replicate');
%         rightresp = imfilter(ROI, rightside, 'replicate');
%         rawresp = nlfilter(ROI,[size(template,1) size(template,2)],fun);
%         rawresp = matchingFiltering(ROI, template, 15);
        rawresp = rawresp.*(1-ROI); %.*rightresp; % favor dark pixels
        respMap(frontBounds(i):end,:,i) = rawresp(1+half_height:end,:);
    end
end
amin = min(respMap(CWL_ROI));
amax = max(respMap(CWL_ROI));
respMap = mat2gray(respMap,[amin amax]);
respMap(~CWL_ROI) = 0;

% visual debug codes
if exist('dispOpt', 'var') && dispOpt
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

