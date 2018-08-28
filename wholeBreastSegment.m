function [segdata] = wholeBreastSegment(dicomdir, savedir, varargin)
%[segdata] = WHOLEBREASTSEGMENT(dicomdir, savedir, varargin)
%   Whole breast segmentation in T1w nonfat-suppressed MRI
%
%   REQUIRED INPUT:
%   - dicomdir: path to the folder from which DICOM files for T1w nonfat-
%     suppressed MRI will be read
%   - savedir: path to the folder in which the segmentation results will be
%     saved
%
%   OPTIONAL INPUT (default values/behaviors will be used if unprovided):
%   - 'Orientation', 'ax'/'sag': imaging orientation of the input MRI; if
%     not provided, the program will try to obtain this information from
%     the DICOM header of the input MRI
%   - 'ID', char string: user assigned ID for the case to process; if
%     undefined, the program with use the 'StudyID' field from the DICOM
%     header
%   - 'v', float/double: half size of the generic template in mm, default
%     value is 2mm, which should work reasonably well in most cases; fine
%     tune can be done around 2mm if needed.
%   - 'q', float/double: 10^q is the weight for the smoothness cost,
%     default value is -1, which should work reasonably well in most cases; 
%     fine tune can be done around -1 if needed.
%   - 'VoxelSize', float/double: isotropic voxel size in mm used during 
%     whole breast segmentation. The recommended values (1mm and 0.5mm)
%     should work well with the default v and q values in most case (the 
%     smaller the voxel size, the longer it takes). If unspecified, the
%     program will determine the voxel size according to the PixelSpacing
%     and SliceThickness info from the DICOM header of the input MRI: if
%     the finest of the two is below 1mm, 0.5mm will be used; otherwise,
%     1mm is used.
%   - 'TemplateDimension', 3 (default) or 2: dimension of the generic and
%     self-adapted templates. 3D templates takes significantly more time 
%     than 2D at a fine voxel resolution (e.g., 0.5mm), but takes almost
%     the same time as 2D at a relative coarse voxel resolution (e.g.,
%     1mm). Meanwhile, 3D templates seem to be more stable and robust than
%     2D. By default 3D templates are used.
%   - 'InterimOutput', true/false(default): if the user wants to save
%     intermediate output, e.g., the ROI for the chest-wall line search, as
%     images, etc. 
%   - 'LateralBounds', true/false: if the user wants to define lateral
%     bounding planes of the breasts, i.e., the two boundaries on the outer
%     sides of the body, to exclude excessive body fat around the armpits.
%     The default is true for sagittal data, and false for axial.
%   - 'BottomBound', true/false(default): if the user wants to define the
%     inferior bounding plane of the breasts to exclude excessive body fat
%     at the belly
%   - 'MedialBounds', true/false(default): if the user wants to define the
%     medial bounding planes of the breasts, i.e., the two boundaries
%     in-between the two breasts.
%   * CAUTION about the breast bounding planes: the breast bounding plane
%     definition algorithms assume excessive non-breast body fat included
%     in the breast MR images. If in reality the field of view of the MRI 
%     is indeed appropriate and the images do not include obviously 
%     excessive non-breast body fat, you should consider avoid defining (at 
%     least some of) the bounding planes. Otherwise you might ending up 
%     with half breasts excluded from the segmented breast mask.
%
%   OUTPUT: segdata: a struct that defines the whole breast segmentation
%           with the following fields (note: this struct actually contains
%           more fields than those listed below; you can safely ignore them
%           for now (either intermediate outputs or reserved for future use))
%     - ID: subject ID
%     - spacing: image voxel resolution: 
%       [yPixelSpacing, xPixelSpacing, sliceThickness]
%     - slices: actual 3D image matrix: [y, x, z]
%     - SeriesDescription: a string extracted from the MRI DICOM header
%       describing the imaging sequence protocol
%     - rightInd: a vector of positive integers (in assending order)
%       corresponding to the part of the input MRI volume that is to the
%       right of the mid-chest. For sagittal data, it contains slice
%       indices; for axial data, it contains column indices.
%     - leftInd: a vector of positive integers (in assending order)
%       corresponding to the part of the input MRI volume that is to the
%       left of the mid-chest. For sagittal data, it contains slice
%       indices; for axial data, it contains column indices.
%     - index1: a vector of positive integers (in assending order)
%       corresponding to the part of the input MRI volume that contains
%       right breast tissue. For sagittal data, it contains slice indices; 
%       for axial data, it contains column indices. index1 can be equal to,
%       or a subset of rightInd.
%     - index2: a vector of positive integers (in assending order)
%       corresponding to the part of the input MRI volume that contains
%       left breast tissue. For sagittal data, it contains slice indices; 
%       for axial data, it contains column indices. index2 can be equal to,
%       or a subset of leftInd.
%     - index (exclusive to axial data): a vector of slice indices (in
%       asscending order) indicating slices that contain breast tissue.
%     - index3 (exclusive to sagittal data): a vector of row indices (in
%       asscending order) indicating rows that contain right breast tissue.
%     - index4 (exclusive to sagittal data): a vector of row indices (in
%       asscending order) indicating rows that contain left breast tissue.
%     - CWL_ROI: a 3D volumetric matrix of type logical and of the same 
%       size as the input MRI volume; true values mark the ROI for the
%       chest-wall line search.
%     - mask1:a 3D volumetric matrix of the same size as the input MRI
%       volume; values above zero comprise the whole breast mask.
%     - VT1/VT2/VT: breast volume (cm^3) for the right, left, and both
%       breasts

% parse optional input
p = inputParser;
addParameter(p, 'Orientation', '', @ischar);
addParameter(p, 'ID', '', @ischar);
addParameter(p, 'v', 2);
addParameter(p, 'q', -1);
addParameter(p, 'VoxelSize', []); % isotropic voxel size in mm
addParameter(p, 'TemplateDimension', 3, @(x)(x==2)||(x==3));
addParameter(p, 'InterimOutput', false, @islogical);
addParameter(p, 'LateralBounds', [], @islogical);
addParameter(p, 'BottomBound', false, @islogical);
addParameter(p, 'MedialBounds', false, @islogical);
parse(p, varargin{:});
Orientation = p.Results.Orientation;
ID = p.Results.ID;
v = p.Results.v;
q = p.Results.q;
VoxelSize = p.Results.VoxelSize;
TemplateDimension = p.Results.TemplateDimension;
InterimOutput = p.Results.InterimOutput;
LateralBounds = p.Results.LateralBounds;
MedialBounds = p.Results.MedialBounds;
BottomBound = p.Results.BottomBound;

% determine image orientation
if ~isempty(Orientation) % if user gives imaging orientation
    if strcmpi(Orientation, 'ax')
        Orientation = 'ax';
    elseif strcmpi(Orientation, 'sag')
        Orientation = 'sag';
    else
        error('Unexpected imaging orientation! The input breast MRI must be either axial or sagittal.'); 
    end
else % try to read out this info from dicom
    dicomlist = dir(fullfile(dicomdir, '*.dcm'));
    if isempty(dicomlist)
        dicomlist = dir(dicomdir);
        dicomlist = dicomlist(3:end);
    end
    % read the first dicom file in the given folder
    dicomfile = dicomlist(1).name;
    dicomheader = dicominfo(fullfile(dicomdir, dicomfile));
    SeriesDescription = dicomheader.SeriesDescription;
    if ~isempty(strfind(SeriesDescription,'ax')) ||  ~isempty(strfind(SeriesDescription,'AX')) || ...
            ~isempty(strfind(SeriesDescription,'Ax')) %#ok<STREMP>
        Orientation = 'ax';
    elseif ~isempty(strfind(SeriesDescription,'sag')) ||  ~isempty(strfind(SeriesDescription,'SAG')) || ...
            ~isempty(strfind(SeriesDescription,'Sag')) %#ok<STREMP>
        Orientation = 'sag';
    else
        error('Failed to identify imaging orientation from DICOM header! Please specify either axial or sagittal orientation.');
    end
end

% call the actual processing function
if strcmpi(Orientation, 'ax')
    if isempty(LateralBounds)
        LateralBounds = false;
    end
    segdata = wholeBreastSegmentAx(dicomdir, savedir, 'v', v, 'q', q, 'VoxelSize', VoxelSize, 'TemplateDimension',...
        TemplateDimension, 'LateralBounds', LateralBounds, 'BottomBound', BottomBound, 'MedialBounds', MedialBounds, ...
        'InterimOutput',InterimOutput, 'ID',ID);
else
    if isempty(LateralBounds)
        LateralBounds = true;
    end
    segdata = wholeBreastSegmentSag(dicomdir, savedir, 'v', v, 'q', q, 'VoxelSize', VoxelSize, 'TemplateDimension',...
        TemplateDimension, 'LateralBounds', LateralBounds, 'BottomBound', BottomBound, 'MedialBounds', MedialBounds, ...
        'InterimOutput',InterimOutput, 'ID',ID);
end

end

