function [validRow] = determineBttmBound_allFOV(masks, slices, spacing, validRow, visualize)
%[validRow] = DETERMINEBTTMBOUND_ALLFOV(masks, slices, spacing, validRow, visualize)
%   此处显示详细说明

% parse input
if ~exist('visualize','var')
    visualize = false;
end

nipple_ht = round(20/spacing(1)); % the empirical lower breast region to be excluded from cut

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

measurement = zeros(size(front_curve,1)-nipple_id-nipple_ht,1);
for i = nipple_id+nipple_ht+1:size(front_curve,1)
    % build a new mask directly connecting from the nipple to the candidate
    % cutting point
    maskx = poly2mask([front_curve(1:nipple_id,1);front_curve(i:end,1);back_curve(:,1)],...
        [front_curve(1:nipple_id,2);front_curve(i:end,2);back_curve(:,2)],...
        size(silhouette,1),size(silhouette,2));
    % get exclusive areas between the two masks
    B = sum( maskx(:) & (~mask(:)) );
    A = sum( mask(:) & (~maskx(:)) );
    measurement(i-nipple_id-nipple_ht) = B - A;
end

% smooth as if a combo of ten 1st-degree polynomials (line segments)
measurement_smooth = smooth(measurement, 10, 'lowess');
if visualize
    figure; plot(measurement,'b'); hold on;
    plot(measurement_smooth,'r'); hold off;
end

[pks,locs,~,p] = findpeaks(max(measurement_smooth)-measurement_smooth, ...
    'NPeaks',2, 'SortStr','descend', 'MinPeakProminence',9);
% P.S.: 'MinPeakProminence' is set to 3^2 in order to be a meaningful
% non-noise minimum area difference
if numel(pks)>1
    pks = pks( p>max(p)/10 );
    locs = locs( p>max(p)/10 );
end
[measurement_sorted,ind] = sort(measurement_smooth);

if any(ismember(measurement_smooth(1:3),measurement_sorted(1:3))) || all(measurement_smooth>0)
    % usally happens in small breasts; then use clustering
    [idx,C]= kmeans([measurement (1:numel(measurement))'],2);
    [~,cmin]= min(C(:,2));
    x = find(idx==cmin, 1, 'last');
elseif isempty(pks)
    [~,x] = min(measurement);
elseif numel(pks)==1
    if ind(1)<=locs
        [~,x] = min(measurement);
    else
        x = locs;
    end
else
    validRow = determineBttmBound_bigFOV(masks, slices, spacing, validRow, visualize);
    return
end

validRow = transpose(front_curve(1:nipple_id+nipple_ht+x,2));

if visualize
    img2show = imoverlayMathWorks(sliceMIP, bwperim(mask), [0 0 1]);
    figure; imshow(img2show); hold on
    % plot(front_curve(:,1),front_curve(:,2), 'r');
    % plot(back_curve(:,1),back_curve(:,2), 'r');
    plot(nipple_xy(1),nipple_xy(2),'yo');
    % the empirical areola region excluded from cut
    plot(front_curve(nipple_id:nipple_id+nipple_ht,1),front_curve(nipple_id:nipple_id+nipple_ht,2),'r');
    % the cutting line
    plot([1 size(sliceMIP,2)],[front_curve(nipple_id+nipple_ht+x,2) front_curve(nipple_id+nipple_ht+x,2)],'g');
    daspect(spacing);
    hold off
end

end

