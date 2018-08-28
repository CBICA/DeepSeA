function [validCol,validSlc] = determineLateralBoundCore(laterality,portion,masks,slices,spacing,validCol,validSlc,visualize)
%[validSlc,validRow]=DETERMINEOUTERBOUNDCORE(laterality,portion,masks,slices,spacing,validCol,validRow[,visualize=false])
%   Determines breast lateral bounds in axial view

% find out whether we'll cut off the end with lower indices, or the other
if strcmpi(laterality, 'right') && strcmpi(portion, 'inner') ...
        || strcmpi(laterality, 'left') && strcmpi(portion, 'outer')
    cutoffLowEnd = false;
elseif strcmpi(laterality, 'right') && strcmpi(portion, 'outer') ...
        || strcmpi(laterality, 'left') && strcmpi(portion, 'inner')
    cutoffLowEnd = true;
else
    error('Invalid laterality or portion indicated. Please check your inputs and try again.');
end

% the empirical nipple region to be excluded from cut; for our data we do
% not find it necessary so we set it to 0
nipple_r = 0; % radius of the nipple region in mm
nipple_r = round(nipple_r/spacing(2)); % convert it to pixel

% find silhouette of breast by mask projection
silhouette = max(masks,[],3);
if cutoffLowEnd % if we are cuting off the end with lower indices
    % flip the silhouette so that we can reuse the codes written for
    % cutting off the other end (i.e., the higher-indexed end)
    silhouette = fliplr(silhouette);
    validCol = size(silhouette, 2) - validCol(end:-1:1) + 1;
end

% find the front surface
front_curve = zeros(numel(validCol), 2);
for j = 1:numel(validCol)
    front_curve(j,:) = [validCol(j) find(silhouette(:,validCol(j)),1,'first')];
end
B = bwtraceboundary(silhouette, fliplr(front_curve(end,:)), 'SW', 8, Inf);
idx = find(B(:,1)==front_curve(1,2) & B(:,2)==front_curve(1,1));
if isempty(idx) || numel(idx)>1
    error('Sth went wrong.');
end
back_curve = fliplr(B(1:idx, :));

% generate a new mask with front_curve and back_curve; this 'mask' is
% different from the 'silhouette' in that it is convex in the front surface
mask = poly2mask([front_curve(:,1);back_curve(:,1)],[front_curve(:,2);back_curve(:,2)],...
    size(silhouette,1),size(silhouette,2));

% identify the nipple as the frontmost point
nipple_y = min(front_curve(:,2));
% in case of multiple frontmost points, use the last one because we are
% cutting off the higher indices
nipple_id = find(front_curve(:,2)==nipple_y, 1, 'last');
nipple_x = front_curve(nipple_id,1);

measurement = zeros(size(front_curve,1)-nipple_id-nipple_r, 1);
for j = nipple_id+nipple_r+1:size(front_curve,1)
    % build a new mask directly connecting from the nipple to the candidate
    % cutting point
    maskx = poly2mask([front_curve(1:nipple_id,1);front_curve(j:end,1);back_curve(:,1)],...
        [front_curve(1:nipple_id,2);front_curve(j:end,2);back_curve(:,2)],...
        size(silhouette,1),size(silhouette,2));
    % get exclusive areas between the two masks
    B = sum( maskx(:) & (~mask(:)) );
    A = sum( mask(:) & (~maskx(:)) );
    measurement(j-nipple_id-nipple_r) = B - A;
    
% 	% codes that generate figure(s) for paper writting
%     im2write = imoverlayMathWorks(im2uint8(mask), maskx(:) & (~mask(:)), [0 1 0]);
%     im2write = imoverlayMathWorks(im2write, mask(:) & (~maskx(:)), [1 0 0]);
%     imwrite(im2write, sprintf('colorArea_%02d.png', j));
%     if j == 115
%         figure; imshow(mat2gray(maskx)); hold on;
%         plot(front_curve(nipple_id+1:j-1, 1), front_curve(nipple_id+1:j-1, 2), 'r--', 'linewidth', 1.25);
%         hold off;
%     end
end

[measurement_sorted,ind] = sort(measurement);
if any(ismember(measurement(1:3),measurement_sorted(1:3))) || all(measurement>0)
    % usally happens in small breasts; then use clustering
    disp('The breast is small and flat; its bounding plane cannot be reliably determined. Operation aborted!');
    if cutoffLowEnd
        % restore original indices (i.e., unflipped) before returning
        validCol = size(silhouette, 2) - validCol(end:-1:1) + 1;
    end
    return;
end

validCol = transpose(front_curve(1:nipple_id+nipple_r+ind(1)-1,1));
criticalCol = validCol(end);

if cutoffLowEnd % restore to unflipped indices, silhouette, etc.
    validCol = size(silhouette, 2) - validCol(end:-1:1) + 1;
    if exist('visualize','var') && visualize
        silhouette = fliplr(silhouette);
        front_curve(:, 1) = size(silhouette, 2) - front_curve(end:-1:1, 1) + 1;
        front_curve(:, 2) = front_curve(end:-1:1, 2);
        nipple_x = size(silhouette, 2) - nipple_x + 1;
        criticalCol = size(silhouette, 2) - criticalCol + 1;
    end
end

% update range of valid rows accordingly as it might be changed as side
% effect
linearInd = find(masks(:,validCol,:));
[~,~,slcSub] = ind2sub(size(masks(:,validCol,:)), linearInd);
validSlc = validSlc( validSlc>=min(slcSub) & validSlc<=max(slcSub) );

if exist('visualize','var') && visualize
%     % codes that generate figure(s) for paper writting
%     figure;
%     imshow(silhouette); hold on;
%     plot(front_curve(:,1),front_curve(:,2), 'r', 'linewidth',1.5);
%     plot(back_curve(:,1),back_curve(:,2), 'g', 'linewidth',1.5);
%     plot(nipple_x,nipple_y,'go', 'markerFaceColor', 'g');
%     daspect(spacing); hold off;    
%     figure;
%     imshow(silhouette); hold on;
%     plot([criticalCol criticalCol], [1 size(slices,1)], 'g--', 'linewidth',1.25);
%     daspect(spacing); hold off;
    
    figure; subplot(1,2,1);
    % mark the critical location in axial view
    imshow(silhouette); hold on;
    plot(front_curve(:,1),front_curve(:,2), 'r');
    plot(nipple_x,nipple_y,'yo');
    plot([criticalCol criticalCol], [1 size(slices,1)],'g'); 
    daspect(spacing); hold off;
    
    % display the first sagittal slice that is EXCLUDED from valid breast
    slices = reformatAx2Sag(slices);
    masks = reformatAx2Sag(masks);
    img2show = imoverlayMathWorks(mat2gray(slices(:,:,criticalCol)), bwperim(masks(:,:,criticalCol)), [0 1 0]);
    subplot(1,2,2);
    imshow(img2show);
    daspect([spacing(3) spacing(1) spacing(2)]);
    
    figure; plot(measurement);
    
%     % codes that generate figure(s) for paper writting
%     figure; plot(front_curve(nipple_id+nipple_r+1:end,1), measurement, 'k', 'linewidth', 1.5); hold on;
%     plot(front_curve(nipple_id+nipple_r+ind(1),1), min(measurement), 'r*', 'MarkerSize', 9); hold off;
end

end

