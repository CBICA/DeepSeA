function []=showLeftbounds(src,event,slices,leftBounds,flagPtr)
%SHOWLEFTBOUNDS 此处显示有关此函数的摘要
%   此处显示详细说明

persistent k; % static to this function
if isempty(k)
    k=1;
end

num_slc = size(slices,3);
if strcmp(event.Key,'escape') % exit signal
    set(flagPtr, 'Value', 1);
    k=1;
    close(src);
    return;
elseif strcmp(event.Key,'leftarrow')
    k = k-1;
    if k<1
        k = 1;
    end
elseif strcmp(event.Key,'rightarrow')
    k = k+1;
    if k>num_slc
        k = num_slc;
    end
elseif strcmp(event.Key,'uparrow')
    k = k-10;
    if k<1
        k = 1;
    end
elseif strcmp(event.Key,'downarrow')
    k = k+10;
    if k>num_slc
        k = num_slc;
    end
elseif strcmp(event.Key,'home')
    k = 1;
elseif strcmp(event.Key,'end')
    k = num_slc;
end

imshow(mat2gray(slices(:,:,k)));
title(['Slice ' num2str(k) '/' num2str(num_slc)]);
if leftBounds(k)>0
    hold on;
    plot([leftBounds(k) leftBounds(k)],[1 size(slices,2)],'g');
    hold off;
end

end

