function [slices,CWL_ROI,leftBounds] = makeISO(slices,spacing,isoSize,CWL_ROI,leftBounds)
%[slices,CWL_ROI,leftBounds]=MAKEISO(slices,spacing,isoSize,CWL_ROI,leftBounds)
%   Detailed explanation goes here

[row,col,slc] = size(slices);

row_iso = round(row*spacing(1)/isoSize);
col_iso = round(col*spacing(2)/isoSize);
slc_iso = round(slc*spacing(3)/isoSize);
[Xq,Yq,Zq] = meshgrid(linspace(1,col,col_iso),linspace(1,row,row_iso),linspace(1,slc,slc_iso));

% interpolate slices: use 'cubic' (more prefered due to better filtering 
% response) or 'linear'; slices produced with 'spline' is visually inferiror
slices = interp3(slices,Xq,Yq,Zq,'cubic');
CWL_ROI = interp3(double(CWL_ROI),Xq,Yq,Zq);
CWL_ROI = CWL_ROI>graythresh(CWL_ROI);
% modify/interpolate frontBounds accordingly
if col_iso ~= col
   leftBounds = round( leftBounds * spacing(2)/isoSize ); 
end
if slc_iso ~= slc
    leftBounds = round( interp1(1:slc, leftBounds, linspace(1,slc,slc_iso)) );
end

end

