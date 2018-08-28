function [CWL_ROI,frontBounds] = processInitMaskAx(mask_init,validSlc,index1,index2)
%[CWL_ROI,frontBounds,validSlc]=PROCESSINITMASK(mask_init)
%   此处显示详细说明

% clear foreground voxels outside the given lateral bounds
toClear = setdiff(1:size(mask_init,2), [index1 index2]);
mask_init(:,toClear,:) = 0;

frontBounds =extractFrontBounds(mask_init,validSlc);
CWL_ROI = defineCWL_ROIwMaskInit(mask_init,frontBounds,validSlc,index1,index2);

end

function [frontBounds] = extractFrontBounds(mask_init, validSlc)
frontBounds = zeros(1,size(mask_init,3));

silhouette = any(mask_init, 3);
% nipple row
tip = find(any(silhouette,2),1);
% find the rough seperation col of two breasts based on half of body voxels
numpxPerCol = sum(sum(mask_init,3),1);
cumNumPxCol = cumsum(numpxPerCol);
cumPrcntPxCol = cumNumPxCol/max(cumNumPxCol);
[~,c] = min(abs(cumPrcntPxCol-.5));
% take front surface of the seperation col as sternum (and landmark)
sternum = find(silhouette(:,c),1);
% [r,c] = find(silhouette);
% butt = min( r(c==min(c)|c==max(c)) );
% [~,butt] = max(sum(silhouette,2));
frontmost = round(tip+(sternum-tip)/3*2);

for k = validSlc
    frontSlcwise = find( any(biggestCC(mask_init(:,:,k)),2), 1 );
    frontBounds(k) = max(frontmost, frontSlcwise);
end

% % uncomment codes below for paper figure generation
% tipx = find(silhouette(tip,:),1);
% rear = find(any(silhouette,2),1,'last');
% rearx = find(silhouette(rear,:),1);
% backend = round((rear+size(silhouette,1))/2);
% figure; imshow(silhouette); hold on;
% plot(tipx,tip, '>', 'markerfacecolor','g','markeredgecolor', 'none','markersize', 8);
% plot(c, sternum, 'v', 'markerfacecolor','g','markeredgecolor','none','markersize',8);
% plot(rearx,rear, '^', 'markerfacecolor','g','markeredgecolor', 'none','markersize', 8);
% plot([1 256], [frontmost frontmost], 'm--','linewidth',1.5);
% plot([1 256], [backend backend], 'm--','linewidth',1.5);
% hold off;

end

function CWL_ROI = defineCWL_ROIwMaskInit(mask_init,frontBounds,validSlc,index1,index2)
% smooth the initial mask a bit
mask_init = imclose(mask_init, ones(3,3,3));

% find out the effective lateral range containing body regions
lateral_lb = min(index1);
lateral_ub = max(index2);

num_validSlc = numel(validSlc);

% extract front/back surface maps; needed in columns where there are no 
% breast, especially in-between two breasts
frontSurfaceMap = zeros(num_validSlc, size(mask_init,2));
backend = 0;
for k = 1:num_validSlc % only slices that contain breast
    slcID = validSlc(k);
    
    % scan through columns to get the front pixel
    for j = lateral_lb:lateral_ub
        i = find(mask_init(:,j,slcID),1,'first');
        if ~isempty(i)
            frontSurfaceMap(k, j) = i;
        end
    end
    
    % find THE pixel on the back
    [r,~] = find(mask_init(:,:,slcID));
    if max(r)>backend
        backend = max(r);
    end

    % slice-wise interpolation of the front surface
%     [~,x,v] = find(frontSurfaceMap(k,:));
%     frontSurfaceMap(k,lateral_lb:lateral_ub) = interp1(x,v, lateral_lb:lateral_ub,'nearest','extrap');
end
% figure; imshow(mat2gray(frontSurfaceMap));
% interpolate / extrapolate where there is no front surface
[rq,cq] = meshgrid(1:num_validSlc,lateral_lb:lateral_ub);
[r, c, v] = find(frontSurfaceMap);
F = scatteredInterpolant(c,r,v,'linear','nearest');
frontSurfaceMap(:,lateral_lb:lateral_ub) = round(F(cq',rq'));
% figure; imshow(mat2gray(frontSurfaceMap));
backend = round( (backend+size(mask_init,1))/2 ); % extend the backend a little bit

CWL_ROI = mask_init;
for k = 1:num_validSlc % only slices that contain breast
    slcID = validSlc(k);
    
    % fill in CWL_ROI between front surface and backend
    for j = lateral_lb:lateral_ub
        CWL_ROI(frontSurfaceMap(k,j):backend,j,slcID) = true;
    end
    
    % remove regions outside the front boundarys for CWL search
    CWL_ROI(1:frontBounds(slcID),:,slcID) = false;    
end

end

