function [mriout] = readmriAx(dicomdir, assignedID)
%READMRIAX Summary of this function goes here
%   Detailed explanation goes here

dicomlist = dir(fullfile(dicomdir, '*.dcm'));
if isempty(dicomlist)
    dicomlist = dir(dicomdir);
    dicomlist = dicomlist(3:end);
end
num_dicom = numel(dicomlist);

fprintf('All %d slices, sorting the %03d', num_dicom, 1);

dicomfile = dicomlist(1).name;
dicomheader = dicominfo(fullfile(dicomdir, dicomfile));
sliceLocations = zeros(num_dicom,1);
sliceLocations(1) = dicomheader.SliceLocation;

for i = 2:num_dicom
    fprintf('\b\b\b%03d', num_dicom, i);
    
    dicomfile = dicomlist(i).name;
    info = dicominfo(fullfile(dicomdir, dicomfile));
    
    % check for data consistency below
    if ~strcmpi(dicomheader.StudyID, info.StudyID)
        error('Dicom files are inconsistent.')
    end
    if ~isequal(dicomheader.PixelSpacing, info.PixelSpacing)
        error('Dicom files are inconsistent.')
    end
    if dicomheader.SliceThickness ~= info.SliceThickness
        error('Dicom files are inconsistent.')
    end
    if ~isequal(dicomheader.ImageOrientationPatient, info.ImageOrientationPatient)
        error('Dicom files are inconsistent.')
    end
    if ~isequal(dicomheader.ImagePositionPatient(1:2), info.ImagePositionPatient(1:2))
        error('Dicom files are inconsistent: ImagePositionPatient.')
    end
    if ~strcmpi(dicomheader.SeriesDescription, info.SeriesDescription)
        error('Dicom files are inconsistent.') 
    end
    if dicomheader.BitDepth ~= info.BitDepth
        error('Dicom files are inconsistent.') 
    end
    
    sliceLocations(i) = info.SliceLocation;
end
fprintf('\n');

if exist('assignedID', 'var')
    studyID = assignedID;
else
    studyID = dicomheader.StudyID;
end

% sort according to SliceLocation
[sliceLocations, I] = sort(sliceLocations, 'ascend');
% check consistency between SliceLocationa and SliceThickness
realthickness = unique(round(diff(sliceLocations),2));
if numel(realthickness)>1
    error('Inconsistent SliceLocation and/or SliceThickness, or missing slice(s).');
end
if abs(realthickness) ~= round(dicomheader.SliceThickness,2)
    error('Inconsistency between SliceLocation and SliceThickness');
end
% allocate memeory for images
slices = zeros(dicomheader.Rows, dicomheader.Columns,num_dicom);
fprintf('All %d slices, reading the %03d', num_dicom, 1);
% read images
for i = 1:num_dicom
    fprintf('\b\b\b%03d', num_dicom, i);
    dicomfile = dicomlist(I(i)).name;
    imgdata = dicomread(fullfile(dicomdir, dicomfile));
    if isempty(imgdata)
        error('DICOM file %s is empty!', dicomfile);
    else
        slices(:,:,i) = imgdata;
    end
end
% in case of left-right reversed image
if dicomheader.ImageOrientationPatient(1)<0
    slices = flip(slices, 2);
    IPP1 = -dicomheader.ImagePositionPatient(1);
else
    IPP1 = dicomheader.ImagePositionPatient(1);
end
% in case of upside down image
if dicomheader.ImageOrientationPatient(5)<0
    slices = flip(slices, 1);
end
fprintf('\n');

% will be used to determine right/left breasts
lrIndicator = double(0:dicomheader.Columns-1)*dicomheader.PixelSpacing(2) + IPP1;

% assign properites
mriout = struct('ID', studyID, ...
    'spacing', [dicomheader.PixelSpacing' dicomheader.SliceThickness], ...
    'slices', slices, ...
    'SeriesDescription', dicomheader.SeriesDescription, ...
    'BitDepth', dicomheader.BitDepth, ...
    'rightInd', find(lrIndicator<0), ...
    'leftInd', find(lrIndicator>=0));

end

