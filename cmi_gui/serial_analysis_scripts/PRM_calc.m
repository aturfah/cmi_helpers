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
