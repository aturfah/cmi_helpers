function data = Load_General(curpwd,C_ID)

% Read file code
str_Col8 = cellfun(@num2str,C_ID(:,8)','un',0).';
num_Col8 = cellfun(@numel,str_Col8','un',0).';    

% Find Framed Data

F_files = sum(cell2mat(num_Col8)==3);
L_files = sum(cell2mat(num_Col8)==6);
R_files = sum(cell2mat(num_Col8)==9);
J_files = sum(cell2mat(num_Col8)==12);


str_R = str_Col8(cell2mat(num_Col8)==9); % list all fixed, label and Reg files from Reg Codes

for j=1:R_files
    str_temp = cell2mat(str_R(j)); % convert cell 2 str
    data.fixed.files(j) = {fullfile(curpwd,C_ID(cell2mat(num_Col8)==3&strcmp(str_Col8,str_temp(end-2:end)),9))}; % Fixed
    data.fixed.label(j) = {fullfile(curpwd,C_ID(cell2mat(num_Col8)==6&strcmp(str_Col8,str_temp(end-5:end)),9))}; % label
    data.moving.files_R(j) = {fullfile(curpwd,C_ID(cell2mat(num_Col8)==9&strcmp(str_Col8,str_temp(1:end)),9))}; % Reg
    data.fixed.mod(j) = C_ID(cell2mat(num_Col8)==3&strcmp(str_Col8,str_temp(end-2:end)),5); % Modality
end
% Identify files for Moving Data and Label
data.moving.files = {fullfile(curpwd,C_ID(cell2mat(num_Col8)==3&~strcmp(str_Col8,str_temp(end-2:end)),9))}; % Moving
data.moving.label = {fullfile(curpwd,C_ID(cell2mat(num_Col8)==6&~strcmp(str_Col8,str_temp(end-5:end)),9))}; % label
data.moving.mod = C_ID(cell2mat(num_Col8)==3&~strcmp(str_Col8,str_temp(end-2:end)),5); % Modality

% Find the Jacobian data sets. Will have 7 or 8 digits
% if J_files>0
%     temp=cell2mat(Jac_Col8);num_temp=numel(temp);
%     data.moving.files_R_Jac = {fullfile(curpwd,C_ID(cell2mat(C_ID(:,8))==str2num(temp(1:num_temp)),9))};
% else
%     data.moving.files_R_Jac = {[]};
% end
% 
%         fprintf(fid, 'Moving Jacobian %5s 2=exists: \n',data.moving.files_R_Jac{:},exist(data.moving.files_R_Jac{:},'file'));
% Load Registered moving label
% if sum(cell2mat(strfind(C_ID(:,7),'label_mR')))>0||sum(cell2mat(strfind(C_ID(:,7),'label_crop_mR')))>0
%     data.moving.files_R_label = {fullfile(curpwd,b{end})};
% else
%     data.moving.files_R_label = {[]};
% end
