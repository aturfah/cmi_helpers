function GTmask_final = Classification_map(data,mask_data,GT_thresholds,cutoff)


maski = {1,2,3,4,5,6}; % Fixed
[d1,d2,d3] = size(data);

% allocate memory
GT = zeros([d1 d2 d3 4]);
datafilt = zeros([d1 d2 d3]);

mask = ismember(mask_data,[maski{:}]); % whole-lung mask
data = data.*mask; % apply whole-lung segmentation to data

h2=waitbar(0,'GT Map: Parallel: 3x3 Median Filter to Exp');movegui(h2,'northwest');waitbar(0,h2);
parfor k = 1:size(data,3);
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