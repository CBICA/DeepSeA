function [diffMap,CWL_ROI,template] = rowwiseLearn(the_row,the_slc,d_field,half_height,half_width,paddedVol,respMap,CWL_ROI,reducedROI,gloTemp)
%[diffMap,CWL_ROI,template]=ROWWISELEARN(the_row,the_slc,d_field,half_height,half_width,paddedVol,respMap,CWL_ROI,reducedROI,gloTemp)
%   Detailed explanation goes here

% IMPORTANT NOTE: it might speed up by padding smooth_slices as a volume in
% the beginning, instead of padding invidual slices inside loops below

w = [0 0 0 0 0];
template = zeros(2*half_height+1,2*half_width+1);
% search for a neareast neighbor in the four upright directions
[nr,ns] = search_4neareastNeighbors(d_field,the_row,the_slc);
for m = 1:4
    if ~isempty(nr{m}) && ~isempty(ns{m})
        % extract template right here
        temp_temp = paddedVol( nr{m}+(0:2*half_height), d_field(nr{m},ns{m})+(0:2*half_width), half_height+ns{m} );
        % weighted by response strength & distance
        w(m) = respMap(nr{m},d_field(nr{m},ns{m}),ns{m}) / (1+abs(the_row-nr{m})+abs(the_slc-ns{m}));
        template = template + w(m)*temp_temp;
    end
end
padded_slice = paddedVol(:,:,half_height+the_slc);
if d_field(the_row,the_slc)
    % if the current location has a valid depth itself, that should
    % be counted as well    
    temp_temp = padded_slice( the_row+(0:2*half_height), d_field(the_row,the_slc)+(0:2*half_width) );
    w(5) = respMap(the_row,d_field(the_row,the_slc),the_slc);
    template = template + w(5)*temp_temp;
    
    % also reduce CWL_ROI to speed up a little bit
    CWL_ROI = CWL_ROI & reducedROI;
end
if any(w)
    template = template/sum(w);
elseif ~isempty(gloTemp)
    warning('CWL_localise:locallyAdpat:noValidNeighbor', ['No valid neighbor ' ...
        'found (including itself) and/or all zeros reponse; using globally adapted template.'])
    template = gloTemp;
else
    error('CWL_localise:locallyAdpat:noValidNeighbor', 'No valid neighbor found (including itself) and/or all zeros reponse.')
end

% row-wise filtering
% first extract effective colums
candi = find(CWL_ROI);
diffMap = ones(size(CWL_ROI));
for j = candi
    jpatch = padded_slice( the_row+(0:2*half_height), j+(0:2*half_width) );
    % calculate difference
    diffMap(j) = norm(template(:)-jpatch(:));
end

end

