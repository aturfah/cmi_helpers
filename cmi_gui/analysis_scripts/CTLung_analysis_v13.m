function output_all = CTLung_analysis_v13(cmiObj,C,GT1orPRM2,HU)
% C = { ID Number(Study Randomization Number) , Patient Initials ,
%       Timepoint(months) , Scan Number , Code , Location }
% HU = {air_ins, trach_ins, blood_ins, air_exp, trach_exp, blood_exp}
% code: 1=moving(exp), 2=fixe(ins), 3=movLab, 4=fixLab, 5=moving Reg,
% 6=Jacobian
% Input:    -cmiObj0: cmi object from Ben's program
%           -C: list of file location
%           -GT1orPRM2: tag to do GT=1, PRM=2, or just calculate T(i)=3
%           -JacTag: 0=no and 1=yes
%           -Need to set PRM Options in CMI program to AllLung (10
%           classifications) as well as cutoff and scatter plot range
%           [-1000 -250].
% Output: -output is a structure with the following fields
% output_all.QAT=GT_thresholds_final; List of QAT thresholds based on Ins
% and Exp
% output_all.ParQAT=output_ParQAT; Values for QAT model [Xexp Xins Yins D]
% output_all.CT=output_final; % Volume (L) and mean HU for Ins and Exp
% output_all.PRM=output_PRM_final; PRM values for -856 and QAT
    % 
% output_all.NRMSE=output_NRMSE; % attemt to determine accuracy of
% registration
% output_all.Dice=output_Dice; % dice coefficient for fixed and registered labels 
% output_all.filename=fname; List of file names

% To copy data:
%   >>prm=cat(1,output_me(:).PRM);prm1=prm(7:7:end,:);  % PRM results
%   >>filename=cat(1,output_me(:).filename);            % filenames
%   >>ct=cat(1,output_me(:).CT);ct1=ct(7:7:end,:);      % CT results
%   >>yInt=cat(1,output_me(:).HUcorr);                  % yIntercept
%   correct
%   open in array viewer


tic;
curpwd = pwd;
N = 3; % this is the column in C that indicates the number of group registrations to perform. Ex different ID or time points.
% try ~ischar(cell2mat(C(:,N)))
%     C_col2 = cell2mat(C(:,N));
%     Cx = unique(C_col2);
%     nID = length(Cx);
% catch
    C_col2=C(:,N);
    nPAT = length(unique(C_col2)); % determine the number of patients in the database
% end

if N==2
    C_col3=num2str(cell2mat(C(:,3)));
    Cx2 = numel((unique(C(:,2))));
    Cx3 = numel((unique(C(:,3))));
    nPAT = length(C_col3)/(Cx2*Cx3);
end
%--------------

% A=1:size(C,1);A=reshape(A,[],nPAT); % generate matrix: number of files x patients eg [5 files 19 patients]

maskj = {1:6}; % Segmetatiofn mask based on lobe segmentation. Will work for fewer segments.
GT_thresholds_final = [];
output_final = [];output_PRM_final = [];output_NRMSE = [];fname={[]};output_Dice = [];output_ParQAT = [];

fixed = [0];
moving = [1];
counti=0;

for i=1:nPAT
    output = [];output_PRM = [];NRMSE=[];Dice_temp = [];ParQAT = [];
    
%         ii=(1+((i-1)*nID):nID*i);
    C_ID=C(:,:);
    % Double check file lists.
        'File Name'
        fname(i,1:2)=[C_ID(1,2),C_ID(1,3)]
        str_ID=[C_ID{1,2},'_',C_ID{1,3}]
        
        if nargin<4||isempty(find(strcmp(HU(:,2),fname(2))==1))
            HUpt = [-1000 -1000 50 -1000 -1000 50];
        else
            'HU Correction Values'
            HUpt = cell2mat(HU(find(strcmp(HU(:,2),fname(2))==1),3:8))
            if isempty(HUpt)
                HUpt = [-1000 -1000 50 -1000 -1000 50]; % Exp then Ins
            end
            HUname = HU(find(strcmp(HU(:,2),fname(2))==1),1:2)
        end
                
    if numel(moving)>1||counti==0
        counti=counti+1;
    end
    
    if sum(cell2mat(C_ID(:,8))==moving(counti))~=0
        % File names:
        a = C_ID(cell2mat(C_ID(:,8))==fixed(counti),end);
        b = C_ID(cell2mat(C_ID(:,8))==moving(counti),end);
        
        %---------------------
        % Copy fixed file names from C_ID
        data.fixed.files = {fullfile(curpwd,a{1})};
        data.fixed.label = {fullfile(curpwd,a{2})};
        fixed_mod = C_ID(cell2mat(C_ID(:,8))==fixed(counti),5);
        data.fixed.mod = fixed_mod(end-1); % remove cells for label
        
       
        % Copy Moving file names from C_ID
        data.moving.files = {fullfile(curpwd,b{1})};
        data.moving.label = {fullfile(curpwd,b{2})};
        data.moving.files_R = {fullfile(curpwd,b{3})};
        
        if strcmpi(data.fixed.mod{1},'exp') %I2E
            data.fixed.HU.values = HUpt(1,4:6); %HU correction for exp
            data.fixed.HU.yInt = -1000-HUpt(1,4);
            data.fixed.HU.slope = 1;
            data.moving.HU.values = HUpt(1,1:3); %HU correction for ins
            data.moving.HU.yInt = -1000-HUpt(1,1);
            data.moving.HU.slope = 1;
        else % E2I only Stanford
            data.fixed.HU.yInt =  -1000-HUpt(1,1); %HU correction for ins
            data.moving.HU.yInt = -1000-HUpt(1,4); %HU correction for exp
        end
        
        % log file names
        %         fprintf(fid, 'Fixed Data %5s 2=exists: %5d \n',data.fixed.files{:},exist(data.fixed.files{:},'file'));
        %         fprintf(fid, 'Fixed Label %5s 2=exists: %5d \n',data.fixed.label{:},exist(data.fixed.label{:},'file'));
        %         fprintf(fid, 'Moving Data %5s 2=exists: %5d \n',data.moving.files{:},exist(data.moving.files{:},'file'));
        %         fprintf(fid, 'Moving Label %5s 2=exists: %5d \n',data.moving.label{:},exist(data.moving.label{:},'file'));
%         fprintf(fid, 'Moving Data_R %5s 2=exists: %5d \n',data.moving.files_R{:},exist(data.moving.files_R{:},'file'));
        
        % Determine if a file is missing
        if ~exist(data.fixed.files{:},'file')||~exist(data.fixed.label{:},'file')||~exist(data.moving.files{:},'file')||...
                ~exist(data.moving.label{:},'file')
            errordlg(['A file was not found for ',C_ID{1,2},'_',C_ID{1,3}]);
            break
        end
        
        % Load Jacobian filename
        if ~isempty(cell2mat(strfind(C_ID(:,7),'Jac')))
            data.moving.files_R_Jac = {fullfile(curpwd,b{4})};
        else
            data.moving.files_R_Jac = {[]};
        end
%         fprintf(fid, 'Moving Jacobian %5s 2=exists: \n',data.moving.files_R_Jac{:},exist(data.moving.files_R_Jac{:},'file'));
        % Load Registered moving label
        if sum(cell2mat(strfind(C_ID(:,7),'label_mR')))>0||sum(cell2mat(strfind(C_ID(:,7),'label_crop_mR')))>0
            data.moving.files_R_label = {fullfile(curpwd,b{end})};
        else
            data.moving.files_R_label = {[]};
        end
        
%         fprintf(fid, 'Moving Label Registerd %5s 2=exists: \n',data.moving.files_R_label{:},exist(data.moving.files_R_label{:},'file'));
%         
%         fclose(fid);
        
        moving_mod = C_ID(cell2mat(C_ID(:,8))==moving(counti),5);
        data.moving.mod = moving_mod(end-2); % remove cells for label and reg
        
        
%         cat(1,data.fixed.files,data.moving.files)
%         cat(1,data.moving.files_R,data.moving.files_R_Jac)
%         cat(1,data.fixed.label,data.moving.label)
%         cat(1,data.fixed.mod,data.moving.mod)
        
        % Save fixed Data to Structure using cmiObj.LoadImg
        for j=1:numel(data.fixed.files)
            cmiObj.loadImg(0,cat(1,data.fixed.files(j),data.fixed.label(j)));
            data.fixed.ct{j,1}=cmiObj.img.mat(:,:,:,1)+data.fixed.HU.yInt; % include HU correction
            data.fixed.label{j,1}=cmiObj.img.mat(:,:,:,2);
            data.fixed.voxel{j,1}=(cmiObj.img.voxsz);
            
            output_Temp=mask_analysis(data.fixed.ct{j,1},data.fixed.label{j,1},...
                data.fixed.voxel{j,1}); % Vol and Mean 4 fixed image
            
            output = cat(2,output,output_Temp);
            local_pwd = pwd;
        end
        % Save moving Data to Structure using cmiObj.LoadImg
        for j=1:numel(data.moving.files)
            cmiObj.loadImg(0,cat(1,data.moving.files(j),data.moving.label(j)));
            data.moving.ct{j,1}=cmiObj.img.mat(:,:,:,1)+data.moving.HU.yInt; % include HU correction
            data.moving.label{j,1}=cmiObj.img.mat(:,:,:,2);
            data.moving.voxel{j,1}=(cmiObj.img.voxsz);
            
            output_Temp=mask_analysis(data.moving.ct{j,1},data.moving.label{j,1},...
                data.moving.voxel{j,1}); % Vol and Mean 4 moving image
            
            output = cat(2,output,output_Temp);
            
            if exist(data.moving.files_R{j},'file')>0
                cmiObj.loadImg(0,cat(1,data.moving.files_R(j)));
                data.moving.reg{j,1}=cmiObj.img.mat(:,:,:,1)+data.fixed.HU.yInt; % include HU correction
            else
                data.moving.reg{j,1}=ones(size(data.fixed.ct{j,1}));
                errordlg('No Registered Data');
            end
                
%             data.moving.reg_label{j,1}=cmiObj.img.mat(:,:,:,2);
            
            % estimate of registration accuracy
            data_gfit=cat(4,data.fixed.ct{1,1},data.moving.reg{1,1});
            mask_gfit=logical(data.fixed.label{1,1});
            A=reshape(data_gfit(repmat(mask_gfit,[1 1 1 2])),[],2);
            NRMSE=cat(2,gfit(A(:,1),A(:,2),'7'),gfit(A(:,1),A(:,2),'4'),goodnessOfFit(A(:,1),A(:,2),'NRMSE'));
            NRMSE
            
            clear data_gfit mask_gfit
            %--------------------------
            
            % Check if Replacement label for Exp (ex. ins_label_R)
            if ~isempty(data.moving.files_R_label{:})
                cmiObj.loadImg(0,cat(1,data.moving.files_R_label(j)));
                data.moving.label_R{j,1}=cmiObj.img.mat;
                % This is where you put DICE coefficient calculation
                A = reshape(data.fixed.label{j,1},size(data.fixed.label{j,1},1),[]);
                B = reshape(data.moving.label_R{j,1},size(data.moving.label_R{j,1},1),[]);
                
                for ii = 1:7
                    if ii<7
                        output_Dice(ii,1) = getDiceCoeff_cjg(ismember(A,ii),ismember(B,ii))
                    else
                        output_Dice(ii,1) = getDiceCoeff_cjg(ismember(A,1:6),ismember(B,1:6))
                    end
                end
            else
                data.moving.label_R{j,1}=[];
            end
            
            % Check if Jac exists. If not set Jac to 1.
            if exist(data.moving.files_R_Jac{j},'file')>0
                cmiObj.loadImg(0,cat(1,data.moving.files_R_Jac(j)));
                data.moving.jac{j,1}=cmiObj.img.mat;
                
                % check if you are using registered label for registered
                % moving image
                if ~isempty(data.moving.files_R_label{:})
                    label = data.moving.label_R{j,1};
                else
                    label = data.fixed.label{j,1};
                end
                
                output_Temp=mask_analysis(data.moving.reg{j,1},label,...
                    data.fixed.voxel{j,1},data.moving.jac{j,1}); % Vol and Mean 4 registered image
                
                output = cat(2,output,output_Temp);
            else
                data.moving.jac{j,1}=ones(size(data.moving.reg{j,1}));
%                 errordlg('No Jacobian Data');
            end
            
            
        end
        cd(curpwd);
        %%
   
        %     p = gcp('nocreate');
        %     if isempty(p)
        %         parpool(1);
        %     end
        
        if GT1orPRM2==1
            str_method = 'Quantitative Air Trapping';
        elseif GT1orPRM2==2
            str_method = 'Parametric Response Mapping';
        else
            str_method = 'Gas Trapping Threshold [T(i)]';
        end
        
        % Check if JacTag is provided
        if nargin<4
            JacTag=0;
        end
        
        % Specify thresholds and cutoffs for Relative Volume (e.g. PRM)
        % calculations
        % PRM Thresholds: This will change with changing GT threshold
        thresh = [2 1 0 -856;1 2 0 -950;1 2 1 -94;1 2 0 -810];
        % cutoff [x y] = [expMR InsF] This will be used for PRM calculations
        cutoff = [1 -1024 -250; 2 -1024 -250];
        % GT Range for Stanford Analysis
        GT_range = cutoff;
        
        h1 = waitbar(0,'Enter the Loop');movegui(h1,'northeast');set(h1,'Name',str_method);
        
        waitbar(i/nPAT,h1,'4/6: Calculate Thresholds');
        % Calculate the GT using T(i) and -856
        for j=1:numel(data.fixed.files)
            if strcmpi(data.fixed.mod{1},'exp') %I2E
                exp = data.fixed.ct{j};
                exp_label = data.fixed.label{j};
                ins = data.moving.ct{j};
                ins_label = data.moving.label{j};
            else % E2I only Stanford
                exp = data.moving.ct{j};
                exp_label = data.moving.label{j};
                ins = data.fixed.ct{j};
                ins_label = data.fixed.label{j};
            end
            
            [T,ParQAT] = Stanford_GT(exp,ins,exp_label,ins_label);
            GT_thresh = cat(2,thresh(1,4),T); % group all thresholds
        end        
        
        if GT1orPRM2==1 % calculate QAT
            waitbar(i/nPAT,h1,'5/6: Calculate QAT');
            % this is where I left off 20160714
            for jj=1:2
                if jj==1 % original exp
                    CT = exp;
                    label = exp_label;
                else % expR and fixed label
                    if strcmpi(data.moving.mod(:),'exp')==1
                        CT = data.moving.reg{1};
                        % determine if you have a replacement Label for
                        % expR
                        if isempty(data.moving.label_R{:})
                            label = data.fixed.label{1};
                        else
                            label = data.moving.label_R{1,1};
                        end
                        jac = data.moving.jac{j,1};
                    else
                        jj=1;
                    end
                end
                x = size(CT,4); % had data where there was 2 Exp for 1 Ins
                
                GT_range_set=GT_range(2,:); % cutoff range is set to [-1024 -250]
                for j=1:x
                    ClassMap = Classification_map(CT,label,GT_thresh(1,1:4),GT_range_set);
                    if jj==1
                        output(:,(jj-1)*4+(7:10)) = Gas_trapping(ClassMap,label,GT_range_set,CT);
                    else
                        output(:,(jj-1)*4+(7:10)) = Gas_trapping(ClassMap,label,GT_range_set,CT,jac);
                    end
                end
            end
            
        elseif GT1orPRM2==2 % calculate PRM
            %% PRM Calculations
            waitbar(i/nPAT,h1,'6/6: Generate PRM Classification Map');
            
            % Fixed Scan
            CTF=data.fixed.ct{strcmpi(data.fixed.mod(:),'Exp'),1};
            if isempty(cell2mat(data.moving.label_R{:}))
                label = data.fixed.label{strcmpi(data.fixed.mod(:),'Exp'),1};
            else
                label = data.moving.label_R{1,1};
            end
            voxel = data.fixed.voxel{strcmpi(data.fixed.mod(:),'Exp'),1};
            jac = ones(size(CTF));
            
            % Registered Scan
            CTR=data.moving.reg{strcmpi(data.moving.mod(:),'Ins'),1};
            
            % TiPharam Bronch study has two Exp. One at T0 and T1. This look
            % sets both Exp as x-axis for PRM using the Ins(T0) as y-axis
            for prmapI = 1:numel(strcmpi(data.moving.mod(:),'Exp'))
                if prmapI == 2
                    CTF=data.moving.reg{strcmpi(data.moving.mod(:),'Exp'),1};
                end
                prmap = PRM_classification(str_ID,local_pwd,cmiObj,label,CTR,CTF,...
                    voxel,GT_thresh(1,1:4),thresh,cutoff);
                x1=size(output,2);
                % output(:,1) is the volume of the fixed image
                output_PRMTemp = PRM_calc(prmap,label,jac,0);
                
                output_PRM = cat(1,output_PRM,output_PRMTemp);
            end
        end
        close(h1);
    else
        GT_thresh = zeros(1,5);
        output = zeros(7,4);
        output_PRM = zeros(7,40);
    end
        %%
        GT_thresholds_final = cat(1,GT_thresholds_final,GT_thresh);
        GT_thresholds = GT_thresh(1,1:4);
        
        output(isnan(output))=1000;output(isinf(output))=0;
        
        output_final = cat(1,output_final,output);
        output_PRM_final = cat(1,output_PRM_final,output_PRM);
        output_NRMSE = cat(1,output_NRMSE,NRMSE);
        output_Dice = cat(1,output_Dice,Dice_temp);
        output_ParQAT = cat(1,output_ParQAT,ParQAT);

               
    ['Elapsed time is ',num2str(round(toc/60,3)),' minutes']
end


% delete(gcp);
cd(curpwd)

output_all.QAT=GT_thresholds_final;
output_all.ParQAT=output_ParQAT;
output_all.CT=output_final;
output_all.PRM=output_PRM_final;
output_all.NRMSE=output_NRMSE;
output_all.Dice=output_Dice; % dice coefficient for fixed and registered labels 
output_all.filename=fname;
output_all.HUcorr=[data.fixed.HU.yInt data.moving.HU.yInt];

assignin('base','output',output_all);
% assignin('base','fname',fname);

['Total elapsed time is ',num2str(round(toc/60,3)),' minutes']

%%
function output=mask_analysis(data,mask_data,voxel,jac)
% calculate the volume and mean HU for each lobe and whole-lung
% Input: exp image, exp label, voxel size
% Output: volume and mean HU

% Specify mask: default is Heidelberg (Right=[1 3], Left=[4 6])
maski = {1:6};
output = zeros(length(maski{1})+1,2);

for j=1:length(maski{1})+1
    % select segmentation
    if j<length(maski{1})+1
        mask_temp = maski{1}(j); % each lobe
    else
        mask_temp = [maski{:}]; % whole-lung
    end
    
    mask = ismember(mask_data,mask_temp);
    if nargin==4
        output(j,1)=sum(jac(mask))*prod(voxel)*1e-6;
    else
        output(j,1)=sum(mask(mask))*prod(voxel)*1e-6;
    end
    output(j,2)=mean(data(mask));
end
clear data mask_data voxel jac mask_temp mask
%%
function [T,ParQAT] = Stanford_GT(exp,ins,mask_exp,mask_ins)

T = zeros(1,3);
maski = {1,2,3,4,5,6}; % Fixed
% Stanford GT: Xexp = 90ile

maskExp = ismember(mask_exp,[maski{:}]);
maskIns = ismember(mask_ins,[maski{:}]);

h2=waitbar(0,'1/2: Stanford QAT: 3x3 Median Filter to Ins');movegui(h2,'northwest');waitbar(0,h2);
% exp = exp.*maskExp;
% ins = ins.*maskIns;

for k = 1:size(maskIns,3);
    waitbar(k/size(maskIns,3),h2,'Filter Ins CT')
    % apply 11x11 median filter to data.*mask(whole-lung)
    insfilt(:,:,k)=medfilt2(ins(:,:,k),[3 3]);
end

waitbar(0,h2,'2/2: Stanford QAT: 3x3 Median Filter to Exp');

for k = 1:size(maskExp,3);
    waitbar(k/size(maskExp,3),h2,'Filter Exp CT')
    % apply 3x3 median filter to data.*mask(whole-lung)
    expfilt(:,:,k)=medfilt2(exp(:,:,k),[3 3]);
end

close(h2)

% identify voxels in mask
expfilt2=expfilt(maskExp);
insfilt2=insfilt(maskIns);

% bypass Goris cutoffs
if nargin<5
    cutoff=[1 -1000 0];
    ['Cutoff: [',num2str(cutoff),']']
end
%--------------

expfilt2(expfilt2<cutoff(1,2)|expfilt2>cutoff(1,3))=[]; % set values outside [output] to zero
insfilt2(insfilt2<cutoff(1,2)|insfilt2>cutoff(1,3))=[]; % set values outside [output] to zero

Xexp = quantile(expfilt2,0.9);
% Stanford GT: X = 90ile and Y = 50ile
[Ins_YX] = quantile(insfilt2,[0.5 0.9]);
Xins = Ins_YX(2);
Yins = Ins_YX(1);
D = Xexp - Xins
for s=1:3
    T(1,s) = Xins-(s-1).*(Xins-Yins)./3-(1-D./343).*(Xins-Yins)/3; % Chest/123/5/May,2003/pg1655
end
ParQAT = [Xexp Xins Yins D];
assignin('base','Stanford',[Xexp,Xins,Yins,D,numel(expfilt2),numel(insfilt2)]);

h = sqrt(Xexp^2+Xins^2);
alpha2 = asind((sind(90)/h)*Xins);
alpha1 = 45-alpha2;
Dnew = sind(alpha1)*(h/sind(90));

[Xexp Xins]
GTnew = -950-(Xins-Xexp);
T(1,4) = GTnew;
clear insF expM mask_fixed mask_moving
%%
function GTmask_final = Classification_map(data,mask_data,GT_thresholds,cutoff)


maski = {1,2,3,4,5,6}; % Fixed
[d1,d2,d3] = size(data);

% allocate memory
GT = zeros([d1 d2 d3 4]);
datafilt = zeros([d1 d2 d3]);

mask = ismember(mask_data,[maski{:}]); % whole-lung mask
data = data.*mask; % apply whole-lung segmentation to data

h2=waitbar(0,'GT Map: Parallel: 3x3 Median Filter to Exp');movegui(h2,'northwest');waitbar(0,h2);
for k = 1:size(data,3);
    % apply 11x11 median filter to data.*mask(whole-lung)
    datafilt(:,:,k)=medfilt2(data(:,:,k),[3 3]);
end

datafilt(datafilt<cutoff(1,2)|datafilt>cutoff(1,3))=0; % set values outside [output] to zero

for jj=1:length(GT_thresholds)
    waitbar(jj/length(GT_thresholds),h2,'Generate mid-filter map.');
    GT_temp = zeros([d1 d2 d3]); % allocate memory
    GT_temp(datafilt<GT_thresholds(1,jj))=1;
    GT(:,:,:,jj) = GT_temp; % final GT masks coded by lobes
end

%Propinquity: apply 15x15 median filter to GTmask (whole-lung)
% GTmaski = {100,200,300,400,500,600}; % Gas trapping by lobe
% GTmask_filt15 = zeros(size(GT));
%
% GTmask = ismember(GT,[GTmaski{:}]); % GTmask for Whole-lung
%
% for kk = 1:size(GTmask,4)
%     GTmask_temp = GTmask(:,:,:,kk);
%     if ~isempty(GTmask_temp(GTmask_temp==1))
%         waitbar(kk/(size(GTmask,4)),h2,'GT Map: Parallel: 15x15 Median Filter to Exp');
%         parfor k = 1:size(GTmask,3)
%             GTmask_filt15(:,:,k,kk) = medfilt2(GTmask(:,:,k,kk),[15 15]);
%         end
%     end
% end
GTmask_filt15 = ones(size(GT));  % used to replace [15 15] filter

GTmask_final = GTmask_filt15.*GT; % keep classification scheme in propinquity map
close(h2)
clear data mask_data GT_thresholds cutoff GTmask_filt15 GT
%%
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
%%
function prmap_final = PRM_classification(str_ID,local_pwd,cmiObj,labelF,ctMR,ctF,voxelF,GT_thresh,thresh,cutoff)

% Propinquity (15x15 median filter on Classification maps) flag
propinquity=0; % 0 = off; 1 = on
if propinquity==1
    str_propinquity = 'Loop through PRM Classifications & Parallel 15x15 Median Filter';
else
    str_propinquity = 'Loop through PRM Classifications. No 15x15 Median Filter';
end

% set cutoff for PRM Options: PRMclass/setOpts=
cmiObj.img.prm.setOpts('cutoff',cutoff);

% Create mask
maski = {1,2,3,4,5,6}; % lobes
mask = ismember(labelF,[maski{:}]); % whole-lung fixed mask

% Pre-Filter MovingR and Fixed
iimg = cat(4,ctF.*mask,ctMR.*mask);
% h2=waitbar(0,'3x3 Wiener Filter for Paired CT Scans');movegui(h2,'northwest');
% for i = 1:size(iimg,4)
%     waitbar(i/2,h2);
%     for k = 1:size(iimg,3)
%         waitbar(k/size(iimg,3),h2,'Wiener Filter Paired CT');
%         % apply 11x11 median filter to data.*mask(whole-lung)
%         iimg(:,:,k,i)=wiener2(iimg(:,:,k,i),[3 3]); % used Matlab Wiener2
%     end
% end

% load MovingR and Fixed: CMIclass/setImg
cmiObj.setImg(iimg,[],size(ctF).*voxelF);

% Load mask: : MaskClass/merge
cmiObj.img.mask.merge('replace',mask); % set cmiObj mask

h2 = waitbar(0,'Loop through GT Thresholds');p=get(h2,'position');
set(h2,'Name','Parametric Response Mapping');p(2)=p(2)-100;
h3=waitbar(0,str_propinquity); set(h3,'position',p);

prmap_final = [];
for j = 1:length(GT_thresh) % loop for calculating PRM per GT threshold
    waitbar(0,h3);
    waitbar(j/(4),h2);
    
    thresh(1,4) = GT_thresh(1,j); % set GT
    thresh(3,4) = thresh(2,4)-thresh(1,4); % recalculate -94HU
    cmiObj.img.prm.setOpts('thresh',thresh); % Run PRM options: PRMclass/setOpts
    
    % Run PRM: ImageClass/calcPRM
    cmiObj.img.calcPRM(2);
    prmap = cmiObj.img.prm.mat;
    prmap_final = cat(4,prmap_final,prmap);
     
    % Topology Analysis (Hoff et al., 2017 Scientific Reports)
    PRMclass = 3; % indicates the PRM classification (Norm=1:2, fSAD=3, emph=4:5, PD=8:10

    if j==1
        mkdir(fullfile(local_pwd,'TopPRM'));
        voxsz=cmiObj.img.voxsz;
        labels=cell2mat(cmiObj.img.labels);
        p = minkowskiFun(prmap,'thresh',PRMclass,'tmode','==','n',[10,10,10],'gridsp',[5,5,5],'voxsz',voxsz,'mask',mask);
        
        scale=[10000 10000 1e5 1e5];
        for MFi=1:4
            img(:,:,:,MFi) = grid2img(p.MF(1,MFi,:),p.ind,p.mask,3,1).*scale(MFi);
        end
        saveMHD(fullfile(local_pwd,'TopPRM',[str_ID,'_',num2str(j),'.mhd']),img,{'V_e4','S_e4','B_e5','X_e5'},size(img(:,:,:,1)).*voxsz);
        saveMHD(fullfile(local_pwd,'TopPRM',[str_ID,'_',num2str(j),'.mhd']),prmap,{'PRM'},size(prmap(:,:,:,1)).*voxsz);
        %---------------------
    end
end

close(h2);close(h3);
clear maski mask labels img prmap_filt
%%
function prm_out = PRM_calc(prmap,labelF,jac,JacTag)

% maski = {1,2,3,4,5,6}; % Gas trapping by lobe
maski = {1:6};
prm_maski = {1,2,3,4,5,6,7,8,9,10};
prm_out1 = zeros(length(maski{1})+1,40);
PRM = zeros(7,4);
maski_Vol = [maski{:}];
num_thresh = size(prmap,4);

% reshape and repmat Jacobian
if JacTag==0
    jac=ones(size(labelF));
end
jac_New = reshape(jac,[],1); % reshape jacobian
jac_New2 = repmat(jac_New,[1,4]); % repmat jacobian vector
clear jac_New

h2 = waitbar(0,'Loop through PRM classifications');movegui(h2,'northwest');
p=get(h2,'position'); p(2)=p(2)-100;set(h2,'Name','Parametric Response Mapping');
h3=waitbar(0,'Loop through Lobes'); set(h3,'position',p);

prm_out1 = [];
% Calculate PRM
for jq = 1:(length(prm_maski)) % loop for PRM class
    waitbar(0,h3);
    waitbar(jq/(length(prm_maski)),h2);
    % reshape and repmat PRM class map
    prm_mask = ismember(prmap,prm_maski{jq}); % acqure prm class mask
    prmap_temp = reshape(prm_mask,[],num_thresh); % reshape and repmat prm_mask
    
    for jqq = 1:(length(maski{1})+1) % loop for lobes
        waitbar(jqq/(length(maski{1})+1),h3);
        if jqq<length(maski{1})+1
            mask = ismember(labelF,maski{1}(jqq));
        else
            mask = ismember(labelF,maski_Vol);
        end
        
        mask4 = repmat(reshape(mask,[],1),[1 num_thresh]); % reshape and repmat lobe mask
        prmap_temp1 = prmap_temp.*mask4; % multiply prm_mask and lobe mask
        
        %         PRM_vol = sum(jac_New2.*prmap_temp1,1).*voxelF.*1e-6;
        %         PRM(jqq,:) = (100.*PRM_vol)./vol(jqq,1);
        
        PRM(jqq,:) = 100.*sum(prmap_temp1,1)./sum(mask4,1);
        
        
    end
    clear prmap_temp1
    prm_out1 = cat(2,prm_out1,PRM);
    
end
prm_out = cat(2,prm_out1(:,1:4:end),prm_out1(:,2:4:end),prm_out1(:,3:4:end),prm_out1(:,4:4:end));
prm_out(isnan(prm_out))=0;prm_out(isinf(prm_out))=0;
close(h2);close(h3);
clear voxel vol jac jac_New2 mask4 prm_mask prmap_temp PRM labelF


