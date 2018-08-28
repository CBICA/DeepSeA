function mriout = readmri( dicomdir, varargin)
% Read the MRI data found in dicomdir. 
% SINTAX:
%     mri = readmri(dicomdir)
%     mri = readmri(dicomdir, 'protocol', 'T1 3D');
%     mri = readmri(dicomdir, 'read_all', true)
%
% INPUT:
%   Required
%     dicomdir <char>   path to the dicom directory
%   Optional
%     protocol <char>   keyword to choose the MR protocol of interest
%                       Default: 'all' for all protocol
%     read_all <bool>   a flag to read in all dicoms in a given path
%                       Default: true
% OUTPUT:
%   mri <struct> data structure in the same standard format as described in
%       the 'Readme.txt'; specifically, the field 'index1' is a vector of 
%       indecies for right-breast slices, and 'index2' left-breast
%
% EXAMPLE:
% >> mri = readmri('.\samples\mri');
% >> slideshow(mri.serie{1});
%
%
%  $Rev: 1633 $:     Revision of last commit
%  $Author: hsiehm@UPHS.PENNHEALTH.PRV $:  Author of last commit
%  $Date: 2018-08-23 10:49:57 -0400 (Thu, 23 Aug 2018) $:    Date of last commit

parseInput

file_list = imlist(dicomdir);
N = length(file_list);
SN = zeros(N, 1);
IN = zeros(N, 1);
SL = zeros(N, 1);
fprintf('Retrieving info:       ')
PN = cell(N, 1);

for n = 1:N
    try
        info = dicominfo(file_list{n});
    catch err
        if strcmp(err.message,'The specified file is not in DICOM format.')
            fprintf('\b\b\b\b\b[%02.0d%%]',floor(100 * n / N))
            continue
        else
            error(err.message);
        end
    end
    if isfield(info, 'SeriesDescription'); PN{n} = info.SeriesDescription;end
    if isfield(info, 'SeriesNumber'); SN(n) = info.SeriesNumber;end
    if isfield(info, 'InstanceNumber'); IN(n) = info.InstanceNumber;end
    if isfield(info, 'SliceLocation'); SL(n) = info.SliceLocation;end
%     if InversionTime ~= info.InversionTime % added by D.Wei on 2016/2/24
%         error('Inverstion time of the designated data is unexpected!')
%     end
    fprintf('\b\b\b\b\b[%02.0d%%]',floor(100 * n / N))
end

% Validating the sequence information
SN = SN(SN~=0);
IN = IN(IN~=0);
PN = PN(~cellfun(@isempty, PN));

% SNu: Unique Series numbers
[SNu,idx] = unique(SN);
PNu = PN(idx);
S = numel(SNu);
fprintf('\n%d files, %d protocols\n', N, S);
fprintf('Retrieving images:     ')
series = repmat(struct('images', [], 'SeriesDescription', []), [S, 1]);
No_images = 0;
remov = false(S, 1);

% S: number of unique series (number)
for s = 1:S
    SN_ref = SNu(s);
    NI_sel = unique(IN(SN == SN_ref));    
    series(s).images = cell(numel(NI_sel), 1);
    N = length(NI_sel);    
    if strcmp(PNu(s), 'localizer')
        remov(s) = true;
        continue
    end
    % N: number of images in series s. Also sort by instance number
    for n = 1:N
        idx = find( (SN == SN_ref) & (IN == NI_sel(n)));
        if numel(idx) > 1;
            error(['More than one slice at instance ', num2str(n)]);
        elseif numel(idx) < 1;
            error(['Unable to find any slice at instance ', num2str(n)]);
        end
        series(s).SeriesDescription = PN{idx};
        series(s).images{n} = file_list{idx};
        series(s).slicelocation(n) = SL(idx);
    end
    No_images = No_images + N;
    fprintf('\b\b\b\b\b[%02.0d%%]',floor(100 * s / S))
end

series(remov) = [];
pflag = (nargin<2);
if ~pflag&&~any(strcmpi(protocol,'all'))
    %Remove images that do not match chosen protocol
    P = length(series);
    remov = false(1, P);
    No_images = 0;
    for p = 1:P        
        remov(p) = ~any(strcmpi(protocol, series(p).SeriesDescription));
        if ~remov(p), No_images = No_images + length(series(p).images);
        end
    end
    series(remov) = [];
    fprintf('\n%d images.\n', No_images)
elseif length(series) == 1;
    fprintf('\n')
else
    %Request user input for a protocol
    fprintf('\n%d images.\n', No_images)
    fprintf('\nLIST OF PROTOCOLS:\n')
    P = length(series);
    for p = 1:P
        fprintf('%d) %s \n', p, series(p).SeriesDescription)
    end
    i = input('Please select protocol: ');
    p = 1:P;
    series(p~=i) = [];
end    

info = dicominfo(series(1).images{1});

mriout = repmat(struct('ID', [], 'spacing', [], 'slices', [], ...
    'mask1',[], 'mask2', [], 'VBD', [], 'FGT',[],'BPE', [],'VT', [],...
    'SeriesDescription', [], 'index1', [], 'index2', [], 'BitDepth', []),size(series));
% 'BitDepth' field added by D.Wei on 2016/2/18 as needed in bias-field correction

for n = 1:length(series)

    %%% Pre and posts dicom images from GE are usually assigned with the same 
    %%% series number and description hence are organized into one directory by 
    %%% sortDicom.m. As we most of the time need just the pre image to be 
    %%% loaded, a GE private tag (0025,1007) and (0025,1011) are used to determine  
    %%% the length of the pre series from the first dicom. In other cases, full 
    %%% length of the series found in a given directory will be used.
    if ~read_all && isfield(info, 'Private_0025_1007') && isfield(info, 'Private_0025_1011');
        hex = flipud(dec2hex(info.Private_0025_1007));
        hex_real = [hex(1,:) hex(2,:) hex(3,:) hex(4,:)];
        numImgInSeries = hex2dec(hex_real);
        hex = flipud(dec2hex(info.Private_0025_1011));
        hex_real = [hex(1,:) hex(2,:)];
        numAcqInSeries = hex2dec(hex_real);
        fprintf('(0025, 1007) Images in Series: %d\n', numImgInSeries);
        fprintf('(0025, 1011) Number of Acquisition: %d\n', numAcqInSeries);
        if numImgInSeries == length(series(n).images);
            numImgInSeries = numImgInSeries/numAcqInSeries;
        elseif numImgInSeries < length(series(n).images); 
            if ~(length(series(n).images)/numImgInSeries == numAcqInSeries);
                numAcqInSeries = length(series(n).images)/numImgInSeries;
                if ~(floor(numAcqInSeries) == ceil(numAcqInSeries));
                    error(['Non-sense value for Number of Acquisition: ' numAcqInSeries]);
                end
                numImgInSeries = length(series(n).images)/numAcqInSeries;
            else
                numImgInSeries = numImgInSeries/numAcqInSeries;
            end
        elseif numImgInSeries > length(series(n).images); 
            warning('Not enough images found for this series in the directory.');
            numImgInSeries = length(series(n).images);
        end
    else
        numImgInSeries = length(series(n).images);
    end
    
    mriout(n).ID = info.AccessionNumber;
    fprintf('Loading images:        ')
    for p = 1:numImgInSeries
        mriout(n).slices(:,:,p) = dicomread(series(n).images{p});
        fprintf('\b\b\b\b\b[%02.0d%%]', floor(100 * p / numImgInSeries))
        if series(n).slicelocation(1)*series(n).slicelocation(p) > 0;
            mriout(n).index1 = [ mriout(n).index1, p ]; % SL < 0 is right breast
        else
            mriout(n).index2 = [ mriout(n).index2, p ]; % SL > 0 is left breast
        end
    end
    fprintf('\n')
    mriout(n).SeriesDescription = series(n).SeriesDescription;
    mriout(n).spacing = [info.PixelSpacing(:); info.SliceThickness];
    mriout(n).BitDepth = info.BitDepth; % 'BitDepth' field added by D.Wei on 2016/2/18 as needed in bias-field correction
end
fprintf('Ok\n')

    function parseInput
        parser=inputParser;
        addRequired(parser, 'dicomdir', @isdir);
        addParameter(parser, 'protocol', 'all', @ischar);
        addParameter(parser, 'read_all', false, @islogical);
        addParameter(parser, 'InversionTime', 1600, @isnumeric);
        
        parse(parser, dicomdir, varargin{:});
        
        dicomdir = parser.Results.dicomdir;
        protocol = parser.Results.protocol;
        read_all = parser.Results.read_all;
        InversionTime = parser.Results.InversionTime;

    end % end of parseInput
    
end % end of readmri

function [image_path] = imlist(folder_path, ext)
% Retrieve list of dicom images in a folder
% SINTAX:
%     image_path = imlist(folder_path)
%     image_path = imlist(folder_path, ext)
% DESCRIPTION:
% Retrieve list of DICOM images in the folder found at
% folder_path. The output, image_path, is a string (or a cell
% array of strings) the full path of the images(s) in that
% folder.
% 
% Said Pertuz
% Nov06/2013

    image_path = [];

    if ~exist(folder_path, 'dir')
        warning('Sorry, I cant find the speficied directory')
        return
    end

    if nargin==2
        flist = dir([folder_path,filesep,'*.',ext]);
        flist([flist.isdir]==1) = [];
        if isempty(flist)
            warning('Sorry, no *.%s files in the speficied directory',ext)
            return
        end    
    elseif nargin==1
        flist = dir(folder_path);
        flist([flist.isdir]==1) = [];
    end

    N = length(flist);
    if N>1, image_path = cell(N, 1);
    else image_path = flist.name; return
    end

    for n = 1:N
        image_path{n} = [folder_path,filesep,flist(n).name];
    end
end
    

