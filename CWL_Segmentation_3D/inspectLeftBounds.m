function inspectLeftBounds(slices,leftBounds)
%INSPECTLEFTBOUNDS(slices,leftBounds)
%   Detailed explanation goes here

flagPtr = libpointer('uint8Ptr',0);
figure('Name',  'Left bounds of CWL_ROIs', ...
    'KeyPressFcn',{@showLeftbounds,slices,leftBounds,flagPtr});
imshow(mat2gray(slices(:,:,1)));
title(['Slice 1/' num2str(size(slices,3))]);
if leftBounds(1)>0 % skip slices that do not contain breast
    hold on;
    plot([leftBounds(1) leftBounds(1)],[1 size(slices,2)],'g');
    hold off;
end
while ~get(flagPtr,'Value')
    pause(1);
end

end

