function prmmap = PRM_classification(cmiObj,dataF,labelF,dataMR,voxelF)

mask = labelF;

% Load Pre-Filter MovingR and Fixed in CMI
iimg = cat(4,dataF.*mask,dataMR.*mask);
cmiObj.setImg(iimg,[],size(dataF).*voxelF);

% Load mask: : MaskClass/merge
cmiObj.img.mask.merge('replace',mask); % set cmiObj mask

% Run PRM: ImageClass/calcPRM
cmiObj.img.calcPRM(2);
prmap = cmiObj.img.prm.mat;
prmmap = {prmap};



