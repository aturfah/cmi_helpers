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