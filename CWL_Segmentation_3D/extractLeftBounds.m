function [leftBounds] = extractLeftBounds(mask_body,validSlc)
%[leftBounds=EXTRACTLEFTBOUNDS(mask_body,validSlc)
%   Detailed explanation goes here

slc = size(mask_body,3);
leftBounds = zeros(1,slc);

silhouette = any(mask_body, 3);
nip = find(any(silhouette,1),1);
[r,c] = find(silhouette);
butt = min( c(r==min(r)|r==max(r)) );
leftmost = round(nip+(butt-nip)/2);

for k = validSlc
    tip = find( any(biggestCC(mask_body(:,:,k)),1), 1 );
    leftBounds(k) = max(leftmost, tip);
end

% % uncomment codes below for paper figure generation
% tipy = find(silhouette(:,nip),1);
% rear = find(any(silhouette,1),1,'last');
% reary = find(silhouette(:,rear),1);
% backend = round((rear+size(silhouette,2))/2);
% figure; imshow(silhouette); hold on;
% plot(nip,tipy, '>', 'markerfacecolor','g','markeredgecolor', 'none','markersize', 16);
% plot(butt, max(r), '>', 'markerfacecolor','g','markeredgecolor','none','markersize',16);
% plot(rear,reary, '<', 'markerfacecolor','g','markeredgecolor', 'none','markersize', 16);
% plot([leftmost leftmost], [1 512], 'm--','linewidth',3);
% plot([backend backend], [1 512], 'm--','linewidth',3);
% hold off;

end

