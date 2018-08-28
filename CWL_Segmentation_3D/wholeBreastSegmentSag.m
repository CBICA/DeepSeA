function [segdata] = wholeBreastSegmentSag(dicomdir,savedir,varargin)
%[segdata]=WHOLEBREASTSEGMENTAX(dicomdir,savedir,varargin)
%   Whole breast segmentation in sagittal T1w nonfat-suppressed MRI; for
%   detailed description of the input and output please see the help info
%   for function "wholeBreastSegment"

visualize = false; % a debug vairable
% Below: preset bottom bound plane detector; the other two options are: 
% (1) determineBttmBound_smlFOV: suitable for breast MRI data that covers
% relatively small range in the superior/inferior direction; and 
% (2) determineBttmBound_bigFOV: suitable for breast MRI data that covers
% relatively big range in the superior/inferior direction.
% The function "determineBttmBound_allFOV" tries to adaptively choose
% either one from the two, but not guaranteed to work correctly all the
% time.
determineBttmBound = @determineBttmBound_allFOV;

% parse input
p = inputParser;
addParameter(p, 'ID', '', @ischar);
addParameter(p, 'v', 2);
addParameter(p, 'q', -1);
addParameter(p, 'VoxelSize', []); % isotropic voxel size in mm
addParameter(p, 'TemplateDimension', 3, @(x)(x==2)||(x==3));
addParameter(p, 'InterimOutput', false, @islogical);
addParameter(p, 'LateralBounds', true, @islogical);
addParameter(p, 'BottomBound', false, @islogical);
addParameter(p, 'MedialBounds', false, @islogical);
parse(p, varargin{:});
ID = p.Results.ID;
v = p.Results.v;
q = p.Results.q;
VoxelSize = p.Results.VoxelSize;
TemplateDimension = p.Results.TemplateDimension;
InterimOutput = p.Results.InterimOutput;
LateralBounds = p.Results.LateralBounds;
MedialBounds = p.Results.MedialBounds;
BottomBound = p.Results.BottomBound;

disp('****************** Begin whole breast segmentation ******************');

% read in the MRI
segdata = readmri(dicomdir);
if ~isempty(ID)
    segdata.ID = ID;
end
segdata.rightInd = segdata.index1;
segdata.leftInd = segdata.index2;

% determine the isotropic voxel size according to the data if not given as
% input
if isempty(VoxelSize)
    finest = min(segdata.spacing);
    if finest < 1
        VoxelSize = .5;
    else
        VoxelSize = 1;
    end
end

disp('Air-breast separation...');
disp('Processing right breast...');
[mask_init1,segdata.index1,segdata.index3] = airBkgSeg(segdata.slices(:,:,segdata.rightInd),segdata.spacing);
disp('Processing left breast...');
[mask_init2,index2,segdata.index4] = airBkgSeg(segdata.slices(:,:,segdata.leftInd),segdata.spacing);
segdata.index2 = index2 + max(segdata.rightInd);
segdata.mask_init = cat(3, mask_init1, mask_init2);
clear mask_init1 mask_init2 index2

if LateralBounds % if requested to determine lateral bouding planes of the breasts
    disp('Determine breast lateral bounds ...');
    disp('Processing right breast...');
    [segdata.index1,segdata.index3] = determineLateralBound('right', 'outer',...
        segdata.mask_init(:,:,segdata.rightInd), segdata.slices(:,:,segdata.rightInd), ...
        segdata.spacing, segdata.index1, segdata.index3, visualize);
    disp('Processing left breast.');
    [segdata.index2,segdata.index4] = determineLateralBound('left', 'outer',...
        segdata.mask_init(:,:,segdata.leftInd), segdata.slices(:,:,segdata.leftInd), ...
        segdata.spacing, segdata.index2-max(segdata.rightInd), segdata.index4, visualize);
    segdata.index2 = segdata.index2 + max(segdata.rightInd);
end

disp('Determine ROI for CWL localization...');
disp('Processing right breast...');
[CWL_ROI1,segdata.leftBounds1] = processInitMask(segdata.mask_init(:,:,segdata.rightInd),...
    segdata.index1,segdata.index3);
disp('Processing left breast.');
[CWL_ROI2,segdata.leftBounds2] = processInitMask(segdata.mask_init(:,:,segdata.leftInd),...
    segdata.index2-max(segdata.rightInd),segdata.index4);
segdata.CWL_ROI = cat(3, CWL_ROI1, CWL_ROI2);
clear CWL_ROI1 CWL_ROI2

if InterimOutput
    % save intermediate results as images for visual debug
    imgdir = fullfile(savedir, segdata.ID, 'CWL_ROI');
    if ~exist(imgdir, 'dir')
        mkdir(imgdir);
    end
    slices = im2uint8(normalize01(segdata.slices));
    for k = 1:size(slices,3)
        sliceRGB = repmat(slices(:,:,k), 1,1,3);
        overlaid = imoverlayMathWorks(sliceRGB, bwperim(segdata.CWL_ROI(:,:,k)), [0,1,0]);
        imwrite( [sliceRGB overlaid], fullfile(imgdir, [num2str(k,'%03.f') '.jpg']) );
    end
    clear imgdir slices sliceRGB overlaid
end

disp('CWL localization...');
disp('Processing 1st (right) breast.');
[CWLs_1,CWL_map1] = CWL_localise(segdata.index1,segdata.index3,segdata.slices(:,:,segdata.rightInd),...
    segdata.spacing,segdata.CWL_ROI(:,:,segdata.rightInd),segdata.leftBounds1,...
    v,q,VoxelSize,TemplateDimension);
% UPDATE CWLs with respect to the range of body mask
CWLs_1 = cutCWL2brst(CWLs_1,segdata.mask_init(:,:,segdata.rightInd),segdata.index1);
% Generate breast mask
segdata.mask1(:,:,segdata.rightInd) = genBrstMask_from_bodyMask_n_CWLmap(...
    segdata.mask_init(:,:,segdata.rightInd),CWL_map1);
disp('Processing 2nd (left) breast.');
segdata.index2 = segdata.index2 - segdata.rightInd(end);
[CWLs_2,CWL_map2] = CWL_localise(segdata.index2,segdata.index4,segdata.slices(:,:,segdata.leftInd),...
    segdata.spacing,segdata.CWL_ROI(:,:,segdata.leftInd),segdata.leftBounds2,...
    v,q,VoxelSize,TemplateDimension);
% UPDATE CWLs with respect to the range of body mask
CWLs_2 = cutCWL2brst(CWLs_2,segdata.mask_init(:,:,segdata.leftInd),segdata.index2);
% Generate breast mask
segdata.mask1(:,:,segdata.leftInd) = genBrstMask_from_bodyMask_n_CWLmap(...
    segdata.mask_init(:,:,segdata.leftInd),CWL_map2);
segdata.index2 = segdata.index2 + segdata.rightInd(end);
% Combine results from both breasts
segdata.CWLs = cat(1, CWLs_1, CWLs_2);
clear CWL_map1 CWLs_1 CWL_map2 CWLs_2

if BottomBound % if requested to determine inferior bouding plane of the breasts
    disp('Determine breast bottom bound...');
    disp('Processing right breast...');
    segdata.index3 = determineBttmBound(segdata.mask1(:,:,segdata.index1),...
        segdata.slices(:,:,segdata.index1),segdata.spacing,segdata.index3,visualize);
    % update segdata fields accordingly
    toClear = setdiff(1:size(segdata.mask1,1), segdata.index3);
    segdata.mask1(toClear,:,segdata.index1) = 0;
    segdata.CWLs(segdata.rightInd) = cutCWL2brst(segdata.CWLs(segdata.rightInd),...
        segdata.mask1(:,:,segdata.rightInd),segdata.index1);
    disp('Processing left breast...');
    segdata.index4 = determineBttmBound(segdata.mask1(:,:,segdata.index2),...
        segdata.slices(:,:,segdata.index2),segdata.spacing,segdata.index4,visualize);
    % update segdata fields accordingly
    toClear = setdiff(1:size(segdata.mask1,1), segdata.index4);
    segdata.mask1(toClear,:,segdata.index2) = 0;
    segdata.CWLs(segdata.leftInd) = cutCWL2brst(segdata.CWLs(segdata.leftInd),...
        segdata.mask1(:,:,segdata.leftInd),segdata.index2-max(segdata.rightInd));
    clear toClear
end

if MedialBounds % if requested to determine medial bouding planes of the breasts
    disp('Determine medial bounds of breasts');
    disp('Processing right breast...');
    [segdata.index1,segdata.index3] = determineLateralBound('right', 'inner',...
        segdata.mask1(:,:,segdata.rightInd), segdata.slices(:,:,segdata.rightInd), ...
        segdata.spacing, segdata.index1, segdata.index3, visualize);
    % update segdata fields accordingly
    toClear = setdiff(segdata.rightInd, segdata.index1);
    segdata.mask1(:,:,toClear) = 0;
    for j = toClear
        segdata.CWLs{j} = [];
    end
    toClear = setdiff(1:size(segdata.mask1,1), segdata.index3);
    segdata.mask1(toClear,:,segdata.index1) = 0;
    segdata.CWLs(segdata.rightInd) = cutCWL2brst(segdata.CWLs(segdata.rightInd),...
        segdata.mask1(:,:,segdata.rightInd),segdata.index1);
    disp('Processing left breast.');
    [segdata.index2,segdata.index4] = determineLateralBound('left', 'inner',...
        segdata.mask1(:,:,segdata.leftInd), segdata.slices(:,:,segdata.leftInd), ...
        segdata.spacing, segdata.index2-max(segdata.rightInd), segdata.index4,visualize);
    segdata.index2 = segdata.index2 + max(segdata.rightInd);
    % update segdata fields accordingly
    toClear = setdiff(segdata.leftInd, segdata.index2);
    segdata.mask1(:,:,toClear) = 0;
    for j = toClear
        segdata.CWLs{j} = [];
    end
    toClear = setdiff(1:size(segdata.mask1,1), segdata.index4);
    segdata.mask1(toClear,:,segdata.index2) = 0;
    segdata.CWLs(segdata.leftInd) = cutCWL2brst(segdata.CWLs(segdata.leftInd),...
        segdata.mask1(:,:,segdata.leftInd),segdata.index2-max(segdata.rightInd));
    clear toClear
end

disp('Quantify breast volume...')
% per-case volume
segdata.VT = calMetrics(segdata.mask1,segdata.spacing);
% per-breast volume
segdata.VT1 = calMetrics(segdata.mask1(:,:,segdata.index1),segdata.spacing);
segdata.VT2 = calMetrics(segdata.mask1(:,:,segdata.index2),segdata.spacing);

disp('Saving results to disk...');
% save results as images for visual debug
imgdir = fullfile(savedir, segdata.ID, 'wholeBreast');
if ~exist(imgdir, 'dir')
    mkdir(imgdir);
end
slices = im2uint8(normalize01(segdata.slices));
for k = 1:size(slices,3)
    sliceRGB = repmat(slices(:,:,k), 1,1,3);
    overlaid = imoverlayMathWorks(sliceRGB, bwperim(segdata.mask1(:,:,k)>0), [0,1,0]);
    imwrite( [sliceRGB overlaid], fullfile(imgdir, [num2str(k,'%03.f') '.jpg']) );
end

% save the segmentation structure to mat file
matdir = fullfile(savedir, segdata.ID);
save(fullfile(matdir, composeSavename(segdata.SeriesDescription, segdata.ID)), 'segdata', '-v7.3');

disp('****************** Whole breast segmentation done! ******************');

end

