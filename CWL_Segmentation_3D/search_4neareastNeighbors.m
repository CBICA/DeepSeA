function [nr,ns] = search_4neareastNeighbors(d_field,cen_r,cen_s)
%[nr,ns] = search_4neareastNeighbors(d_field,cen_r,cen_s)
%   此处显示详细说明

% initiate
nr = cell(1,4);
ns = cell(1,4);

%% search in the same row as the current central location
d_row = d_field(cen_r, :);
valid_slc = find(d_row);
slc_diff = valid_slc - cen_s;

m_diff = min(slc_diff(slc_diff>0));
if ~isempty(m_diff)
    ns{1} = cen_s + m_diff;
    nr{1} = cen_r;
end

m_diff = max(slc_diff(slc_diff<0));
if ~isempty(m_diff)
    ns{2} = cen_s + m_diff;
    nr{2} = cen_r;
end

%% search in the same slice as the current central location
d_slc = d_field(:, cen_s);
valid_row = find(d_slc);
row_diff = valid_row - cen_r;

m_diff = min(row_diff(row_diff>0));
if ~isempty(m_diff)
    nr{3} = cen_r + m_diff;
    ns{3} = cen_s;
end

m_diff = max(row_diff(row_diff<0));
if ~isempty(m_diff)
    nr{4} = cen_r + m_diff;
    ns{4} = cen_s;
end

end

