function Output = serial_analysis_v2(cmiObj,C)
% C = { ID Number(Study Randomization Number) , Patient Initials ,
%       Timepoint(months) , Scan Number , Code , Location }
% HU = {air_ins, trach_ins, blood_ins, air_exp, trach_exp, blood_exp}
% code: 1=moving(exp), 2=fixe(ins), 3=movLab, 4=fixLab, 5=moving Reg,
% 6=Jacobian
% Input:    -cmiObj: cmi object from Ben's program
%           -C: list of file location
%           -GT1orPRM2: tag to do GT=1, PRM=2, or just calculate T(i)=3
%           -JacTag: 0=no and 1=yes
%           -Need to set PRM Options in CMI program to AllLung (10
%           classifications) as well as cutoff and scatter plot range
%           [-1000 -250].
% Output: -output is a cell with the following fields
%   Output_Col2: -name of 2nd column in database
%   Output_Col3: -name of 3rd column in database
%   Output_PRM:  -left over and PRM values
%       Ex. 4 color PRM = [0 1 2 3 4 5 6 7]; 0 is voxels not classified
%   Output_tPRM: -scaled mean tPRM values over mask [Ve4 Se4 Be5 Xe5]

% To copy data:
%   >>data = cell2mat(output(:,3))
%   open in array viewer

tic;

C(:,1:3) = cellfun(@num2str,C(:,1:3)','un',0).'; % convert all cells to strings
nPAT = [length(C(:,1)) length(unique(C(:,1))) length(unique(C(:,2))) length(unique(C(:,3)))]; % determine the number of unique folders


N = reshape(1:nPAT(1),nPAT(1)/prod(nPAT(2:4)),[]);

fixed = [0];
moving = [1];

curpwd = pwd;
Output = {[]};

for i = 1:size(N,2)
    
    %% Reading Files from Database
    C_ID = C(N(:,i),:);
    data = Load_General(curpwd,C_ID);
    
    %% Loading Data from Files in Database Using CMI Program
    % Save fixed Data to Structure using cmiObj.LoadImg
    for j=1:numel(data.fixed.files)        
        cmiObj.loadImg(0,cat(1,data.fixed.files{j},data.fixed.label{j}));
        data.fixed.ct{j,1}=cmiObj.img.mat(:,:,:,1); % no HU correction
        
        data.fixed.label{j,1}=cmiObj.img.mat(:,:,:,2);
        data.fixed.voxel{j,1}=(cmiObj.img.voxsz);

        local_pwd = pwd;
    end
    
    % Save moving Data (1. Moving & Label, 2. Moving_R, 3. Jac, 4.Label_R) to Structure using cmiObj.LoadImg
    
    for j=1:numel(data.moving.files)
        
        % 1. Moving Data and Label
        cmiObj.loadImg(0,cat(1,data.moving.files{j},data.moving.label{j}));
        data.moving.ct{j,1}=cmiObj.img.mat(:,:,:,1); % no HU correction
        data.moving.label{j,1}=cmiObj.img.mat(:,:,:,2);
        data.moving.voxel{j,1}=(cmiObj.img.voxsz);
        
        % 2. Registered Moving data
        if exist(data.moving.files_R{j}{1},'file')>0
            cmiObj.loadImg(0,cat(1,data.moving.files_R{j}));
            data.moving.reg{j,1}=cmiObj.img.mat(:,:,:,1); % no HU correction
        else
            data.moving.reg{j,1}=ones(size(data.fixed.ct{j,1}));
            errordlg('No Registered Data');
        end
        
        % 3. Check if Jac exists. If not set Jac to 1.
        if isfield(data.moving,'files_R_Jac')
            cmiObj.loadImg(0,cat(1,data.moving.files_R_Jac(j)));
            data.moving.jac{j,1}=cmiObj.img.mat;
        else
            data.moving.jac{j,1}=ones(size(data.moving.reg{j,1}));
            %                 errordlg('No Jacobian Data');
        end
        
        % 4. Check if Replacement label for Exp (ex. ins_label_R)
        if isfield(data.moving,'files_R_label')
            cmiObj.loadImg(0,cat(1,data.moving.files_R_label(j)));
            data.moving.label_R{j,1}=cmiObj.img.mat;
        else
            data.moving.label_R{j,1}={[]};
        end
        
    end
    cd(curpwd);
    
    %% Functions for PRM & tPRM
    % set range (-1024 -250) and thresholds (4 colors) for PRM Options: PRMclass/setOpts=
    cmiObj.img.prm.setOpts('cutoff',[1 -1024 -250; 1 -1024 -250]);
    cmiObj.img.prm.setOpts('thresh',[2 1 0 -856; 1 2 0 -950; 1 2 1 -94]);
    prm_range = 0:cmiObj.img.prm.nprm;
    
    % Create Fixed mask used in PRM and tPRM
    maski = {1,2,3,4,5,6}; % lobes
    mask = ismember(data.fixed.label{:},[maski{:}]); % whole-lung fixed mask
    
    %---Calculate PRM
    data.moving.PRM = PRM_classification(cmiObj,data.fixed.ct{:},mask,data.moving.reg{:},...
        data.fixed.voxel{:});
    
    %---Calculate tPRM
%     data.moving.tPRM_Maps = tPRM(data.moving.PRM{:},3,mask,data.fixed.voxel{:});
    
    %% Outputs
    %---Calculate Mask Volumes and Mean Values for Fixed and Moving
    Output_Fixed = mask_analysis(data.fixed.ct{:},data.fixed.label{:},data.fixed.voxel{:});
    Output_Moving = mask_analysis(data.moving.ct{:},data.moving.label{:},data.moving.voxel{:});
    
    %---Calculate PRM Percent Volumes
    prm = data.moving.PRM{:};
    a=prm(prm>=0);
    Output_PRM = 100*histc(a,prm_range)'./sum(histc(a,prm_range));
    
    %---Calculate tPRM over Mask
%     Output_tPRM = mean(reshape(data.moving.tPRM_Maps(repmat(mask,[1 1 1 size(data.moving.tPRM_Maps,4)])),[], size(data.moving.tPRM_Maps,4)),1);
    Output_tPRM = 0;
    
    Output_Values = cat(2,Output_Fixed,Output_Moving,Output_PRM,Output_tPRM);
    Output(i,1:3) = cat(2,C_ID(1,2),C_ID(1,3),Output_Values)
end
assignin('base','Output',Output);
assignin('base','Output_Values',Output_Values);
openvar('Output_Values')

toc;
