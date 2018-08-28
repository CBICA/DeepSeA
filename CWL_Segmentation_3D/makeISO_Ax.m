function [slices,CWL_ROI,frontBounds] = makeISO_Ax(slices,spacing,isoSize,CWL_ROI,frontBounds)
%[slices,CWL_ROI,leftBounds]=MAKEISO_AX(slices,spacing,isoSize,CWL_ROI,leftBounds)
%   Assuming pixel is square within image plane

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
if row_iso ~= row
    frontBounds = round( frontBounds * spacing(1)/isoSize );
end
if slc_iso ~= slc
    frontBounds = round( interp1(1:slc, frontBounds, linspace(1,slc,slc_iso)) );
end

end

