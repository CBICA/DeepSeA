function d = dMap(slices,spacing,CWL_ROI,leftBounds,filter_size_half,lambdaAdjuster,iso_spacing,Tdimension,visualize)
%d=DMAP(slices,spacing,CWL_ROI,leftBounds,filter_size_half,lambdaAdjuster,iso_spacing,visualize)
%   This function assumes inputs are only from effective slices selected by
%   the user

[row,col,slc] = size(slices);

if ~exist('iso_spacing', 'var') || isempty(iso_spacing)
    iso_spacing = min(spacing);
end
% make isotropic if necessary
if any(spacing~=iso_spacing)
    fprintf('Interpolating for isotropic voxels ... ');
    [slices,CWL_ROI,leftBounds] = makeISO(slices,spacing,iso_spacing,CWL_ROI,leftBounds);
    disp('Done!');
end
[row_iso,col_iso,slc_iso] = size(slices);

% visual debug codes
if exist('visualize', 'var') && visualize
    inspectLeftBounds(slices,leftBounds);
end

%% initial response and corresponding CWL
% convert size of filter from physical dimension to pixels
filter_size_half = round(filter_size_half/iso_spacing);
if filter_size_half < 2
    filter_size_half = 2; % at minimum there must be 2 pixels to construct the filter
end
fprintf('Generically filtering ... ');
if Tdimension == 3
    [respMap,smooth_slices] = CWL_likelihood_filtering3D(slices,CWL_ROI,leftBounds,filter_size_half,visualize);
else
    [respMap,smooth_slices] = CWL_likelihood_filtering(slices,CWL_ROI,leftBounds,filter_size_half,visualize);
end
disp('Done!');

% extract CWL map from response map
if exist('visualize','var') && visualize
    CWL_map = extractCWLmap_byMaxResp(respMap, slices);
else
    CWL_map = extractCWLmap_byMaxResp(respMap);
end

% % thresh = graythresh(respMap(CWL_ROI));
% % strongMap = respMap<thresh;

%% local template learning and difference map generation
disp('Self adapting:');
if Tdimension == 3
    globalAdaptFcn = @globalAdaptSag3d;
    localAdaptFcn = @rowwiseLearn3d;
else
    globalAdaptFcn = @globalAdaptSag;
    localAdaptFcn = @rowwiseLearn;
end

half_width = 2*filter_size_half; % width of the template patch
half_height = 2*filter_size_half; % height of the template patch
paddedVol = padarray(smooth_slices, [half_height half_width half_height], 'replicate');

fprintf('Globally adapting ... ');
gloTemp = globalAdaptFcn(CWL_map,paddedVol,respMap,half_height,half_width);
% % codes that generate figure for paper
% imwrite(gloTemp, 'global_template.png');
disp('done!');

% produce localised per-row template and then filter each row
fprintf('Locally adapting ... ');
[d_field] = CWL_mapTo_d_field(CWL_map);
reducedROI = imdilate(CWL_map, true(1,1+2*half_width)); % later used to reduce ROI to speed up a little bit
diffMap = ones(row_iso,col_iso,slc_iso);
parfor k = 1:slc_iso % parfor can be used here by design
    for i = 1:row_iso
        [diffMap(i,:,k),CWL_ROI(i,:,k)]= localAdaptFcn(i,k,d_field,half_height,half_width,...
            paddedVol,respMap,CWL_ROI(i,:,k),reducedROI(i,:,k),gloTemp);        
%         % codes that generate figure for paper
%         if k==55 && i==133
%             imwrite(template, 'local_template.png');
%         end
    end
end
amin = min(diffMap(CWL_ROI));
amax = max(diffMap(CWL_ROI));
diffMap = mat2gray(diffMap,[amin amax]);
diffMap(~CWL_ROI) = 1;
disp('done!');
% visually debug codes
if exist('visualize', 'var') && visualize
    implay(1-diffMap);
end
disp('Self adaptation completed!');

%% extract another set of d0, optimize it and return result
% initialize the depth map for energy minimization
[~,d0] = min(diffMap,[],2);
% extract true positives as the biggest connected component
CWL_map = false(row_iso,col_iso,slc_iso);
[slcSub,rowSub] = meshgrid(1:slc_iso,1:row_iso);
CWL_map(sub2ind([row_iso,col_iso,slc_iso],rowSub(:),d0(:),slcSub(:)))= true;    
d0 = CWL_mapTo_d_field(biggestCC(CWL_map));
% initiate false positive locations with random true positives
ind_true = find(d0);
ind_doubt = find(~d0);
d0(~d0) = d0( ind_true(randi(numel(ind_true),1,numel(ind_doubt))) );

disp('Begin depth field optimization ... ');
d = dMapEnergy(d0,diffMap,lambdaAdjuster,visualize);

% visually debug codes
if exist('visualize', 'var') && visualize
    CWL_map = false(row_iso,col_iso,slc_iso);
    [slcSub,rowSub] = meshgrid(1:slc_iso,1:row_iso);
    CWL_map(sub2ind([row_iso,col_iso,slc_iso],rowSub(:),round(d(:)),slcSub(:)))...
        = true;
    slices_4disp = mat2gray(slices);
    slices_4disp(CWL_map) = 1;
    implay(slices_4disp);
end

% restore original dimensions if necessary
if slc_iso ~= slc || row_iso ~= row
    [X,Y] = meshgrid( linspace(1,slc,slc_iso), linspace(1,row,row_iso) );
    [Xq,Yq] = meshgrid(1:slc,1:row);
    d = interp2(X,Y, d, Xq,Yq, 'linear');
end
if col_iso ~= col
    d = d*col/col_iso;
end

end

