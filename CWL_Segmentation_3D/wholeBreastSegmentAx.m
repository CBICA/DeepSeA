function [segdata] = wholeBreastSegmentAx(dicomdir,savedir,varargin)
%[segdata]=WHOLEBREASTSEGMENTAX(dicomdir,savedir,varargin)
%   Whole breast segmentation in axial T1w nonfat-suppressed MRI; for
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
addParameter(p, 'LateralBounds', false, @islogical);
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
segdata = readmriAx(dicomdir);
if ~isempty(ID)
    segdata.ID = ID;
end

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
[segdata.mask_init,segdata.index,segdata.index1,segdata.index2]=airBkgSegAx...
    (segdata.slices,segdata.rightInd,segdata.leftInd,segdata.spacing);

if LateralBounds % if requested to determine lateral bouding planes of the breasts
    disp('Determine breast outer bounds...');
    [segdata.index,segdata.index1,segdata.index2] = determineLateralBoundAx('outer',...
        segdata.mask_init,segdata.slices,segdata.spacing,segdata.rightInd,segdata.leftInd,...
        segdata.index,segdata.index1,segdata.index2, visualize);
end

disp('Determine ROI for CWL localization...');
[segdata.CWL_ROI,segdata.frontBounds] = processInitMaskAx(...
    segdata.mask_init,segdata.index,segdata.index1,segdata.index2);

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
% filter_size = 2*v+1mm, lambda = 10^q
[segdata.CWLs,CWL_map] = CWL_localiseAx(segdata.index,segdata.index1,segdata.index2,...
    segdata.slices,segdata.spacing,segdata.CWL_ROI,segdata.frontBounds,...
    v,q,VoxelSize,TemplateDimension); 
segdata.mask1 = genBrstMask_from_bodyMask_n_CWLmapAx(segdata.mask_init,...
    CWL_map,segdata.index1,segdata.index2);
clear CWL_map

if BottomBound % if requested to determine inferior bouding plane of the breasts
    disp('Determine breast bottom bound...');
    segdata.index = determineBttmBoundAx(segdata.mask1(:,[segdata.index1 segdata.index2],:),...
        segdata.slices(:,[segdata.index1 segdata.index2],:),segdata.spacing,...
        segdata.index,determineBttmBound, visualize);
    % update segdata fields accordingly
    % clear masks and CWLs below bottom bound
    toClear = setdiff(1:size(segdata.mask1,3),segdata.index);
    segdata.mask1(:,:,toClear) = 0;
    for i = toClear
        segdata.CWLs{i} = [];
    end
end

if MedialBounds % if requested to determine medial bouding planes of the breasts
    disp('Determine medial bounds of breasts...');
    [segdata.index,segdata.index1,segdata.index2] = determineLateralBoundAx('inner',...
        segdata.mask_init,segdata.slices,segdata.spacing,segdata.rightInd,segdata.leftInd,...
        segdata.index,segdata.index1,segdata.index2, visualize);
    % update segdata fields accordingly
    % clear masks and CWLs below bottom bound
    toClear = setdiff(1:size(segdata.mask1,3),segdata.index);
    segdata.mask1(:,:,toClear) = 0;
    for i = toClear
        segdata.CWLs{i} = [];
    end
    % clear between-breasts region in mask1
    toClear = setdiff(1:size(segdata.mask1,2), [segdata.index1 segdata.index2]);
    segdata.mask1(:,toClear,:) = 0;
end

disp('Quantify breast volume...')
voxelVolume = prod(segdata.spacing/10);
segdata.VT = voxelVolume * sum(segdata.mask1(:)>0);
segdata.VT1 = voxelVolume *sum(sum(sum( segdata.mask1(:,segdata.index1,segdata.index)>0 )));
segdata.VT2 = voxelVolume *sum(sum(sum( segdata.mask1(:,segdata.index2,segdata.index)>0 )));

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

