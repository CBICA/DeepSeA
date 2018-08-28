function [validRow] = determineBttmBound_bigFOV(masks, slices, spacing, validRow, visualize)
%[validRow] = DETERMINEBTTMBOUND_BIGFOV(masks, slices, spacing, validRow, visualize)
%   此处显示详细说明

% find silhouette of breast by mask projection
masks = masks>0;
silhouette = max(masks,[],3);
% maxium projection slice
sliceMIP = normalize01(max(slices,[],3));

% find breast front surface on the silhouette
front_curve = zeros(numel(validRow),2);
for i = 1:numel(validRow)
    c = find(silhouette(validRow(i),:),1,'first');
    if ~isempty(c) % sometimes the given validRow does not match the real mask
        front_curve(i,:) = [c validRow(i)];
    end
end
front_curve(front_curve(:,1)==0,:) = []; % remove empty rows
B = bwtraceboundary(silhouette, fliplr(front_curve(end,:)), 'NE', 8, Inf, 'counterclockwise');
idx = find(B(:,1)==front_curve(1,2) & B(:,2)==front_curve(1,1));
back_curve = fliplr(B(1:idx,:));

% generate a new mask with front_curve and back_curve; this 'mask' is
% different from the 'silhouette' in that it is convex in the front surface
mask = poly2mask([front_curve(:,1);back_curve(:,1)],[front_curve(:,2);back_curve(:,2)],...
    size(silhouette,1),size(silhouette,2));

% locate nipple
[~,nipple_id]= min(front_curve(:,1));
nipple_xy = front_curve(nipple_id,:);

% devide front_curve into 3 segments by clustering
[idx,C] = kmeans(front_curve,3);
[~,idy]= sort(C(:,2));
% get the mid seg which must be breast
front_curve_brst = front_curve(idx==idy(2),:);

% now look at the lower portion
front_curve_lower = front_curve(idx==idy(3),:);
[idx,C] = kmeans(front_curve_lower,2);
[~,idy]= min(C(:,2));
% the roi within which to look for the cut point
front_curve_roi = front_curve_lower(idx==idy,:);

measurement = zeros(size(front_curve_roi,1),1);
for i = 1:size(measurement)
    % locate current point in the entire front curve
    j = find(front_curve(:,2)==front_curve_roi(i,2));
    if isempty(j) || numel(j)>1
        error('Sth wrong with front_curve_roi.')
    end
    % build a new mask directly connecting from the nipple to the candidate
    % cutting point
    maskx = poly2mask([front_curve(1:nipple_id,1);front_curve(j:end,1);back_curve(:,1)],...
        [front_curve(1:nipple_id,2);front_curve(j:end,2);back_curve(:,2)],...
        size(silhouette,1),size(silhouette,2));
    % get exclusive areas between the two masks
    B = sum( maskx(:) & (~mask(:)) );
    A = sum( mask(:) & (~maskx(:)) );
    measurement(i) = B - A;
end

% necessry to smoothing a little?
% measurement = smooth(measurement, 25, 'rloess');

% [measurement_sorted,ind]= sort(measurement, 'descend');
% if any(ismember(measurement(1:3),measurement_sorted(1:3)))
%     % may happen in relative big breasts when pressed hard
%     [~,locs] = findpeaks(smooth(measurement, round(numel(measurement)/smoothSegments), 'rloess'), ...
%         'NPeaks',2,'SortStr','ascend');
%     xmax = locs(1);
% else
%     xmax = ind(1);
% end
[measurement_sorted,idx] = sort(measurement);
if any(ismember(measurement(1:3),measurement_sorted(1:3))) || all(measurement>0)
    % usally happens in small/flat breasts; then use clustering
    [idx,C]= kmeans([measurement (1:numel(measurement))'],2);
    [~,idy]= min(C(:,2));
    x = find(idx==idy, 1, 'last');
else
    x = idx(1);
end

boundRow = front_curve_roi(x,2);
validRow = front_curve(1,2):boundRow;

if exist('visualize','var') && visualize
    figure; plot(measurement);
    img2show = imoverlayMathWorks(sliceMIP, bwperim(mask), [0 0 1]);
    figure; imshow(img2show); hold on
    % plot(front_curve(:,1),front_curve(:,2), 'r');
    % plot(back_curve(:,1),back_curve(:,2), 'r');
    plot(nipple_xy(1),nipple_xy(2),'yo'); % nipple
    plot(front_curve_brst(:,1),front_curve_brst(:,2),'r');
    plot(front_curve_roi(:,1),front_curve_roi(:,2),'y');
    % the cutting line
    plot([1 size(sliceMIP,2)],[boundRow boundRow],'g');
    daspect(spacing);
    hold off
end

end

