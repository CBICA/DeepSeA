function [CWL_ROI,leftBounds,validSlc,validRow] = defineCWL_ROI_wLeftbounds(mask_init,leftBounds,validSlc,validRow)
%[CWL_ROI]=DEFINECWL_ROI_WLEFTBOUNDS(mask_init,leftBounds,validSlc,validRow)
%   Detailed explanation goes here

% parse input
if ~exist('validRow','var')
    [r,~] = find(any(mask_init(:,:,validSlc),3));
    validRow = unique(r);
end

linearInd = find(mask_init);
[~,colSub,~] = ind2sub(size(mask_init), linearInd);
backend = round( (max(colSub)+size(mask_init,2))/2 );

CWL_ROI = false(size(mask_init));
for k = validSlc % only slices that contain breast
    [r,c]=find(mask_init(:,:,k));
    if ~any(c>=leftBounds(k))
        % The breast region is totally to the left of the left boundary of
        % CWL ROI. This can happen in <extremely> large breasts in which
        % an iosolated region of breast appear in front of the body but the
        % undetected CWL region arealdy fades/darkens. In this case we
        % exclude this slice from breast region in accordance with our
        % clinician generated ground truth.
        leftBounds(k) = 0;
        validSlc = setxor(validSlc, k);
        continue
    end
    
    CWL_ROI(validRow,leftBounds(k):backend,k) = true;
    % remove air outside (to the left of) the front breast surface
    r=r(c>=leftBounds(k));
    c=c(c>=leftBounds(k)); % now only need to consider the portion to the right of the left boundary
    ur = unique(r);
    % the first effective row in the initial mask and above
    frontSentinel = min(c(r==ur(1)));
    CWL_ROI(validRow(1):ur(1),leftBounds(k):frontSentinel-1,k) = false;
    % all rows in the middle
    for i = 2:numel(ur)-1
        frontSentinel = min(c(r==ur(i)));
        CWL_ROI(ur(i),leftBounds(k):frontSentinel-1,k) = false;
    end
    % the last effective row in the initial mask and below
    frontSentinel = min(c(r==ur(end)));
    CWL_ROI(ur(end):validRow(end),leftBounds(k):frontSentinel-1,k) = false;
end

% update the validRow
[r,~] = find(any(mask_init(:,:,validSlc),3));
validRow = transpose(unique(r));

end

