function [validSlc,validCol1,validCol2] = determineLateralBoundAx(portion,masks,slices,spacing,rightInd,leftInd,validSlc,validCol1,validCol2,visualize)
%[validSlc,validCol1,validCol2]=DETERMINELATERALBOUNDAX(portion,masks,slices,spacing,rightInd,leftInd,validSlc,validCol1,validCol2,visualize)
%   INPUT
%   portion: 'outer' / 'inner'

% parse inputs
if ~exist('visualize', 'var')
    visualize = false;
end

% two breasts are processed seperately

% right breast first
[validCol1,validSlc1] = determineLateralBoundCore('right',portion,masks(:,rightInd,:),...
    slices(:,rightInd,:),spacing,validCol1,validSlc,visualize);

% then left breast
[validCol2,validSlc2] = determineLateralBoundCore('left',portion,masks(:,leftInd,:),...
    slices(:,leftInd,:),spacing,validCol2-max(rightInd),validSlc,visualize);
validCol2 = validCol2 + max(rightInd);

% union valid slices range
validSlc = union(validSlc1, validSlc2);

end

