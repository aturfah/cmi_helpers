function output=mask_analysis(data,mask_data,voxel)
% claculate the volume and mean HU for whole-lung
% Input: exp image, exp label, voxel size
% Output: volume and mean HU

% Specify mask: default is Heidelberg (Right=[1 3], Left=[4 6])
maski = {1:6};
mask = ismember(mask_data,[maski{:}]);

output(1,1)=sum(mask(mask))*prod(voxel)*1e-6;
output(1,2)=mean(data(mask));
