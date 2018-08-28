function inspectFrontBounds(slices,frontBounds)
%INSPECTFRONTBOUNDS(slices,frontBounds)
%   Detailed explanation goes here

flagPtr = libpointer('uint8Ptr',0);
figure('Name',  'Front bounds of CWL_ROIs', ...
    'KeyPressFcn',{@showFrontbounds,slices,frontBounds,flagPtr});
imshow(mat2gray(slices(:,:,1)));
title(['Slice 1/' num2str(size(slices,3))]);
if frontBounds(1)>0 % skip slices that do not contain breast
    hold on;
    plot([1 size(slices,2)],[frontBounds(1) frontBounds(1)],'g');
    hold off;
end
while ~get(flagPtr,'Value')
    pause(1);
end

end

