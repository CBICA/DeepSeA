function [CWL_ROI,leftBounds,validSlc,validRow] = processInitMask(mask_init,validSlc,validRow)
%[ABBs,mask_body,CWL_ROI,leftBounds]=PROCESSINITMASK(mask_init,validSlc,validRow)
%   此处显示详细说明

% % obtain mask_body
% ABBs = cell(size(mask_init,3),1);
% mask_body = mask_init;
% for k = validSlc
%     ABBs{k} = extractABBfromBrstMask(mask_init(:,:,k));
%     mask_body(:,:,k) = genMask1fromABB(ABBs{k},size(mask_init,1),size(mask_init,2));
% end

leftBounds =extractLeftBounds(mask_init,validSlc);
[CWL_ROI,leftBounds,validSlc,validRow] = defineCWL_ROI_wLeftbounds(mask_init,leftBounds,validSlc,validRow);
% mask_body = mask_init | CWL_ROI;
% CWL_ROI = mask_body;
% for k = validSlc
%     CWL_ROI(:,1:leftBounds(k),k) = false;
% end

% % Remove armpit ghost by examing air-breast interface
% quarterRow = round(size(mask_init,1)/4);
% for k = validSlc
%     mask_snpt = mask_body(1:quarterRow,:,k);
%     air_bound_x = zeros(quarterRow,1);
%     for i = 1:quarterRow
%         air_bound_x(i)=find(mask_snpt(i,:),1);
%     end
%     
%     frontSentinel=max(air_bound_x);
%     r=find(air_bound_x==frontSentinel,1,'last');
%     if r>1
%         CWL_ROI(1:r-1,1:frontSentinel,k) = false;
%         mask_body(1:r-1,1:frontSentinel,k) = false;
%     end
%     % TO DO: better front boundary extrapolation after armpit detection
%     % - Tried 3D extrapolation with interpolants and non-parametric surface 
%     %   fitting, results were bad
% end

% upperVolMask = mask_body(1:quarterRow,:,validSlc); % the use of validSlc 
% % insures every included slice contains breast
% frontSurfDepthMap = mask2frontSurfMap(upperVolMask);
% [~,axBound] = max(flipud(frontSurfDepthMap),[],1);
% axBound = quarterRow+1-axBound;
% [X,Y] = meshgrid(1:numel(validSlc),1:quarterRow);
% axMap = Y< repmat(axBound,quarterRow,1);
% if any(axMap(:)) % it is only necessary when armpit ghost is detected
%     x = X(~axMap);
%     y = Y(~axMap);
%     sf = fit([x,y], frontSurfDepthMap(~axMap), 'lowess');
%     frontSurfDepthMap(axMap) = round( sf(X(axMap),Y(axMap)) );
% %     frontSurfDepthMap(axMap)=griddata(x,y,frontSurfDepthMap(~axMap),X(axMap),Y(axMap),'v4');
%     frontSurfDepthMap(frontSurfDepthMap<1) = 1;
%     frontSurfDepthMap(frontSurfDepthMap>col) = col;
% end
% for k = validSlc
%     for i = 1:quarterRow
%         mask_body(i,1:frontSurfDepthMap(i,k-min(validSlc)+1),k) = false;
%         CWL_ROI(i,1:frontSurfDepthMap(i,k-min(validSlc)+1),k) = false;
%     end
% end

end

