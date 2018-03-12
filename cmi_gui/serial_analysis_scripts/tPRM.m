function img = tPRM(prmap,PRMclass,mask,voxelF)
 % Topology Analysis (Hoff et al., 2017 Scientific Reports)

p = minkowskiFun(prmap,'thresh',PRMclass,'tmode','==','n',[10,10,10],'gridsp',[5,5,5],'voxsz',voxelF,'mask',logical(mask));

scale=[10000 10000 1e5 1e5];
for MFi=1:4
    img(:,:,:,MFi) = grid2img(p.MF(1,MFi,:),p.ind,p.mask,3,1).*scale(MFi);
end
% saveMHD(fullfile(local_pwd,'TopPRM',[str_ID,'_',num2str(j),'.mhd']),img,{'V_e4','S_e4','B_e5','X_e5'},size(img(:,:,:,1)).*voxsz);
% saveMHD(fullfile(local_pwd,'TopPRM',[str_ID,'_',num2str(j),'.mhd']),prmap,{'PRM'},size(prmap(:,:,:,1)).*voxsz);
