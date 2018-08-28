function d = dMapAx(slices,spacing,CWL_ROI,frontBounds,filter_size_half,lambdaAdjuster,iso_spacing,Tdimension,visualize)
%d=DMAPAX(slices,spacing,CWL_ROI,frontBounds,filter_size_half,lambdaAdjuster,iso_spacing,visualize)
%   This function assumes inputs are only from effective slices selected by
%   the user

[row,col,slc] = size(slices);

if ~exist('iso_spacing', 'var') || isempty(iso_spacing)
    iso_spacing = min(spacing);
end
% make isotropic if necessary
if any(spacing~=iso_spacing)
    fprintf('Interpolating for isotropic voxels ... ');
    [slices,CWL_ROI,frontBounds] = makeISO_Ax(slices,spacing,iso_spacing,CWL_ROI,frontBounds);
    disp('Done!');
end
[row_iso,col_iso,slc_iso] = size(slices);

% visual debug codes
if exist('visualize', 'var') && visualize
    inspectFrontBounds(slices,frontBounds);
end

%% initial response and corresponding CWL
% convert size of filter from physical dimension to pixels
filter_size_half = round(filter_size_half/iso_spacing);
if filter_size_half < 2
    filter_size_half = 2; % at least there must be 2 pixels to construct the filter
end
fprintf('Generically filtering ... ');
if Tdimension == 3
    [respMap,smooth_slices] = CWL_likelihood_filteringAx3D(slices,CWL_ROI,frontBounds,filter_size_half,visualize);
else
    [respMap,smooth_slices] = CWL_likelihood_filteringAx(slices,CWL_ROI,frontBounds,filter_size_half,visualize);
end
disp('Done!');

% extract CWL map from response map
if exist('visualize','var') && visualize
    CWL_map = extractCWLmap_byMaxRespAx(respMap, slices);
else
    CWL_map = extractCWLmap_byMaxRespAx(respMap);
end

% % thresh = graythresh(respMap(CWL_ROI));
% % strongMap = respMap<thresh;

%% local template learning and difference map generation
disp('Self adapting:');
if Tdimension == 3
    globalAdaptFcn = @globalAdaptAx3d;
    localAdaptFcn = @colwiseLearn3d;
else
    globalAdaptFcn = @globalAdaptAx;
    localAdaptFcn = @colwiseLearn;
end

half_width = 2*filter_size_half; % width of the template patch
half_height = 2*filter_size_half; % height of the template patch
paddedVol = padarray(smooth_slices, [half_height half_width half_width], 'replicate');

fprintf('Globally adapting ... ');
gloTemp = globalAdaptFcn(CWL_map,paddedVol,respMap,half_height,half_width);
disp('done!');

% produce localised per-column template and then filter each column
fprintf('Locally adapting ... ');
D0 = CWL_mapToDepthFieldAx(CWL_map);
reducedROI = imdilate(CWL_map, true(1+2*half_height,1)); % will be used to reduce ROI to speed up a little bit
diffMap = ones(row_iso,col_iso,slc_iso);
parfor k = 1:slc_iso % parfor can be used here by design
    for j = 1:col_iso
        [diffMap(:,j,k),CWL_ROI(:,j,k)]= localAdaptFcn(j,k,D0,half_height,half_width,...
            paddedVol,respMap,CWL_ROI(:,j,k),reducedROI(:,j,k),gloTemp);
    end
end
amin = min(diffMap(CWL_ROI));
amax = max(diffMap(CWL_ROI));
diffMap = mat2gray(diffMap,[amin amax]);
diffMap(~CWL_ROI) = 1;
disp('done!');
% visually debug codes
if exist('visualize','var') && visualize
    implay(1-diffMap);
end
disp('Self adaptation completed!');

%% extract D^1, optimize it and return result
% initialize the depth map for energy minimization
[~,D1] = min(diffMap,[],1);
% d1 = reshape(d1,col_iso,slc_iso); % reshape d0 into a depth field
% if any(d1(:)==1)
%     warning('Found row(s) of all ones in the difference map!')
%     d1(d1==1) = row_iso;
% end
CWL_map = false(row_iso,col_iso,slc_iso);
[slcSub,colSub] = meshgrid(1:slc_iso,1:col_iso);
CWL_map(sub2ind([row_iso,col_iso,slc_iso],D1(:),colSub(:),slcSub(:)))...
    = true;
D1 = CWL_mapToDepthFieldAx(biggestCC(CWL_map));
% initiate false positive locations with random true positives
ind_true = find(D1);
ind_doubt = find(~D1);
D1(~D1) = D1( ind_true(randi(numel(ind_true),1,numel(ind_doubt))) );
% % initiate false positive locations by interpolation/extropolation
% [y,x] = find(d0);
% F = scatteredInterpolant(x,y,d0(d0>0), 'natural');
% [yq,xq] = find(~d0);
% d0(~d0) = F(xq,yq);

disp('Begin depth field optimization ... ');
d = dMapEnergyAx(D1, diffMap, lambdaAdjuster, visualize);

% visually debug codes
if exist('visualize','var') && visualize
    CWL_map = false(row_iso,col_iso,slc_iso);
    [slcSub,colSub] = meshgrid(1:slc_iso,1:col_iso);
    CWL_map(sub2ind([row_iso,col_iso,slc_iso],round(d(:)),colSub(:),slcSub(:)))...
        = true;
    slices_4disp = mat2gray(slices);
    slices_4disp(CWL_map) = 1;
    implay(slices_4disp);
end

% restore original dimensions if necessary
if slc_iso~=slc || col_iso~=col
    [X,Y] = meshgrid( linspace(1,slc,slc_iso), linspace(1,col,col_iso) );
    [Xq,Yq] = meshgrid(1:slc,1:col);
    d = interp2(X,Y, d, Xq,Yq, 'linear');
end
if row_iso ~= row
    d = d*row/row_iso;
end

end

