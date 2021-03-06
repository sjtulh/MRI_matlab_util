
% Simultanes Phase Unwrapping and Removal of chemical Shift (SPURS) via Graph Cut

% [wwater wfat wfreq wunwph_uf unwphw N_std] = spurs_gc(iField,TE,CF,voxel_size,SUBSAMPLE)
%% HOW TO USE
% SPURS with full-matrix: [wwater wfat wfreq wunwph_uf unwphw N_std] = spurs_gc(iField,TE,CF,voxel_size);
% SPURS with subsample: [wwater wfat wfreq wunwph_uf unwphw N_std] = spurs_gc(iField,TE,CF,voxel_size,2);
% Note that if the water fat map totally swap, try conj(iField) instead of iField as input
%% output: 
%     - wwater: the water map 
%     - wfat:   the fat map 
%     - wfreq:  the field map in rad after running IDEAL as fine tunning,input for QSM
%     - wunwph_uf:  the field map after unwrapping and unfat,initial guess for IDEAL
%     - unwphw: phase unwrapping result

%% input: 
%      - iField : a multi-echo 4 dimentional data (Note that if the water
%      fat map totally swap, try conj(iField) instead of iField as input)
%      - how to choose iField or conj(iField) as input: 
%             if PrecessionIsClockwise = 1, [] = spurs_gc(conj(iField),TE,CF,voxel_size);
%             if PrecessionIsClockwise = -1, [] = purcs_gc(iField,TE,CF,voxel_size);

% written by Jianwu Dong  2014.2.10
% last modified by Jianwu 2014.9.3
% last modified by Jianwu, add voxel_size when computing the gradient
%%
%   Modify initial guess of field map for IDEAL:
%       Subtract k*2pi in order to make it as close to 0 as possible
% ==== By Zhe Liu, 3/16/2017 ==== %

function [wwater wfat wfreq wunwph_uf unwphw N_std ] = spurs_gc_mod(iField,TE,CF,voxel_size,SUBSAMPLE,dfat,Mask)
if nargin<5
    SUBSAMPLE = 1;
end

energy = [];
iField0 = iField;
[sx sy sz necho] = size(iField);

if abs((TE(2)-TE(1))-(TE(3)-TE(2)))< 0.0002
    [iFreq_raw N_std] = Fit_ppm_complex(iField);
else
    [iFreq_raw N_std] = Fit_ppm_complex_TE(iField,TE);
end

iFreq_raw(isnan(iFreq_raw))=0;
iFreq_raw(isinf(iFreq_raw))=0;
iFreq_raw0 = iFreq_raw(:,:,:);
iFreq_raw1 = iFreq_raw(:,:,:);

if SUBSAMPLE == 2
    iField = iField(1:SUBSAMPLE:end,1:SUBSAMPLE:end,:,:);
    iFreq_raw1 = iFreq_raw0(1:SUBSAMPLE:end,1:SUBSAMPLE:end,:);
    N_std1 = N_std(1:SUBSAMPLE:end,1:SUBSAMPLE:end,:);
    Mask = Mask(1:SUBSAMPLE:end,1:SUBSAMPLE:end,:);
    voxel_size(1) = voxel_size(1)/0.5;
    voxel_size(2) = voxel_size(2)/0.5;
end

iMag = sqrt(sum(abs(iField).^2,4));



delta_TE = TE(2) - TE(1);

if nargin < 6 || isempty(dfat)
    dfat = -3.5e-6*CF;
end
    dyna_range = 1/delta_TE;
    effect_fat_Hz = dfat + floor( (0.5*dyna_range-dfat)/dyna_range)*dyna_range;
    effect_fat_rad = effect_fat_Hz/dyna_range*2*pi;


p = 2;
w1 = effect_fat_rad/pi


if (w1 > 0)
    w = w1;
    [unwphw,iter,erglist] = phase_unwrap_3d(iFreq_raw1,p,iMag,voxel_size);
    energy = erglist;
    % Correct 2*pi shift
    unwphw = correct_2pi_shift(unwphw, Mask);
    [wkappa,wm_fat,wunwph_uf,iter,erglist,wkiter] = unwrap_unfat_3dP(voxel_size,iMag,w,unwphw,p); % a small mistake found here. 2014.2.21
    energy = [energy erglist];
end


if (w1 < 0)
    w = -w1;  
    [unwphw,iter,erglist] = phase_unwrap_3d(iFreq_raw1,p,iMag,voxel_size);
    energy = erglist;
    % Correct 2*pi shift
    unwphw = correct_2pi_shift(unwphw, Mask);
    [wkappa,wm_fat,wunwph_uf,iter,erglist,wkiter] = unwrap_unfat_3dN(voxel_size,iMag,w,unwphw,p);
    energy = [energy erglist];
end

% Correct 2*pi shift
[wunwph_uf, n_2pi] = correct_2pi_shift(wunwph_uf, Mask);
wkappa = wkappa - n_2pi*2;

% interpolation to get the field map
if SUBSAMPLE == 2
    allX = 1:sx;
    allY = 1:sy;
    subX = 1:SUBSAMPLE:sx;
    subY = 1:SUBSAMPLE:sy;
    [ALLX,ALLY] = meshgrid(allY(:),allX(:));
    [SUBX,SUBY] = meshgrid(subY(:),subX(:));
    fm = zeros(sx,sy,sz);
    for ind = 1:sz
        fm(:,:,ind) = interp2(SUBX,SUBY,wunwph_uf(:,:,ind),ALLX,ALLY,'*spline');
    end 
    
    k_vals = unique(wkappa);
    k_app = (fm - iFreq_raw0)/pi;

    for i = 1:length(k_vals)
        dk(:,:,:,i) = abs(k_app - k_vals(i));
    end 
    [a b] = min(dk,[],4);

    K = b;
    for i = 1:length(k_vals)
        K(b==i)= k_vals(i);
    end
    wunwph_uf = iFreq_raw + K*pi;
 end

% iField = iField0;
iField = conj(iField0);     % Fix bug, to be consistent with Fit_ppm_complex
% 

%% Run IDEAL take wunwph_uf as initial guess 
%% 

% Modify wunuph_uf: subtract multiple of 2*pi
% wunwph_uf_raw = wunwph_uf;
% tmp = wunwph_uf;
% idxs_try = [-5:5];
% clear tmps
% for i = 1:length(idxs_try)
%     tmps(i) = abs(sum(col((tmp - 2*pi*idxs_try(i)).*Mask)));
% end
% [~,idx_tmp] = min(tmps);
% idx_wrap = idxs_try(idx_tmp)
% wunwph_uf = tmp - 2*pi*idx_wrap;



if 1
    [wwater wfat wfreq] = fit_IDEAL(single(iField(:,:,:,:)), TE, dfat, (wunwph_uf)/(2*pi*delta_TE),[],5);
    if sum(abs(wfat(:)).^2)>sum(abs(wwater(:)).^2)
        disp(['potential water fat swap']);
        wunwph_uf = wunwph_uf+effect_fat_rad;
        [wwater wfat wfreq] = fit_IDEAL((iField(:,:,:,:)), TE, dfat, (wunwph_uf)/(2*pi*delta_TE),[],5);
    end    
else
    [xx yy zz] = size(wunwph_uf);
    R2s = zeros([1 xx*yy*zz]);
    % note that when include R2s in IDEAL, the result may contain may noisey point
    [wwater wfat wfreq R2s] = fit_IDEAL_R2((iField(:,:,:,:)), TE, dfat, (wunwph_uf - 4*pi)/(2*pi*delta_TE),R2s,2);
    
    %% maybe try the following: but need choosing the filter parameter of hann_low for different dataset
%   [wwater wfat wfreq R2s] = fit_IDEAL_R2(conj(iField(:,:,:,:)), TE, dfat, (wunwph_uf)/(2*pi*delta_TE),R2s,5);
%    R2s(R2s>100)=0;
%    r2sf = hann_low(R2s, voxel_size, 80);
%    r2sf = real(r2sf);
%    [wwater wfat wfreq r2s] = fit_IDEAL_R2(conj(iField(:,:,:,:)), TE, dfat, (wunwph_uf)/(2*pi*delta_TE),r2sf,0);   
end

wfreq = wfreq*2*pi*delta_TE;
