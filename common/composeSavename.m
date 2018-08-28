function [savename] = composeSavename(seriesDescription,pid)
%[savename]=COMPOSESAVENAME(seriesDescription,pid)
%   Detailed explanation goes here

seriesDescription(~isstrprop(seriesDescription, 'alphanum')) = '_';
savename = sprintf('%s_%s', pid, seriesDescription);

end

