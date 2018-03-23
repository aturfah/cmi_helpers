function GT = Gas_trapping(GT_map,mask_data,cutoff,exp,jac)

maski = {1:6};
GT = zeros(length(maski{1})+1,4);

mask = ismember(mask_data,[maski{:}]);
'GT Range Used for QAT'
cutoff

for i=1:4
    GT_map1 = GT_map(:,:,:,i);
    GT_map1(exp(mask)<cutoff(1,2)|exp(mask)>cutoff(1,3)) = 0; % set voxels outside Exp[-985 -250] to 0
    GT_map2(:,:,:,i) = GT_map1;
end

if nargin==5
    str = 'Jac Corrected';
else
    str = 'ExpM';
end

h2 = waitbar(0,[str, ': Calculating the GT = Vol of GT/ Vol Lobe']);movegui(h2,'northwest');
% Calculate GT
for jq = 1:(length(maski{1})+1)
    waitbar(jq/(length(maski)+1),h2);
    if jq<(length(maski{1})+1)
        mask = ismember(mask_data,maski{1}(jq));
    else
        mask = ismember(mask_data,[maski{:}]);
    end
    
    GT_map_temp = reshape(GT_map2,[],4);
    
    mask_temp = repmat(reshape(mask,[],1),[1,4]);
    
    if nargin==5
        jac(exp(mask)<cutoff(1,2)|exp(mask)>cutoff(1,3)) = 0; % set voxels outside cutoff to 0
        jac_New = repmat(reshape(jac,[],1),[1,4]);
        GT_vol = sum(jac_New.*GT_map_temp.*mask_temp,1);
        vol = sum(reshape(jac,[],1).*reshape(mask,[],1)); % sum of all jac in mask and cutoff.
    else
        GT_vol = sum(GT_map_temp.*mask_temp,1);
        vol = numel(exp(exp(mask)>cutoff(1,2)&exp(mask)<cutoff(1,3))); % numel between cutoff.
    end
    
    GT(jq,1:4) = (100.*GT_vol)./vol;
end
close(h2);
clear GTmask_final voxel vol jac