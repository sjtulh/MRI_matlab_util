% Prepare preconditioner for Total Field Inversion (TFI)
%   [P] = gen_precond(Mask, Ps, R2star, R2star_th, R2star_P)
%
%   output
%   P - Precondioner:   Used as element-wise multiplier
%   
%   input
%   Mask - a binary 3D matrix denoting the Region Of soft tissue
%          typically obtained by thresholding the magnitude
%   Ps (optional) - weight for strong source in a binary preconditioner
%                   (given 1 for weak source)
%                   Default: 30
%   R2star (optional) - use the R2* for extracting "strong" soft tissue
%   R2star_th (optional) - threshold level for R2*,
%                          ( Mask of "strong" soft tissue := R2* >= R2star_th )
%                           Default: 30 sec^-1
%   R2star_P (optional) - weight for high R2*
%
%
%   When using the code, please cite 
%     1. Liu, Z., et al. (2016). "Preconditioned total field inversion (TFI) 
%           method for quantitative susceptibility mapping." Magnetic Resonance 
%           in Medicine DOI: 10.1002/mrm.26331.
%     2. Liu, Z., et al. (2017). "Optimization of Preconditioned Total Field 
%           Inversion for Whole head QSM and Cardiac QSM." ISMRM 2017.
%           p3663
%
% ==== by Zhe Liu, 5/8/2017 ==== %

function [P, Mask_G] = gen_precond(Mask, Ps, R2star, R2star_th, R2star_P)

if nargin < 2
    Ps = 30;
end

if nargin < 5
    R2star_P = Ps;
end

if nargin < 4
    R2star_th = 30;
end

if nargin < 3
    flag_R2s = 0;
else
    flag_R2s = 1;
end

if flag_R2s
    % Use R2*
    Mask_weak = (R2star < R2star_th) & Mask;
    Mask_strong = ~Mask_weak;
    Mask_R2s = R2star >= R2star_th;
    P = 1*single(Mask_weak) + Ps*single(~Mask) + R2star_P*single(Mask & ~Mask_weak);
    Mask_G = 1*single(Mask) + 1/Ps*single(Mask);
else
    % Do not use R2*
    P = 1*single(Mask) + Ps*single(~Mask);
    Mask_G = 1*single(Mask) + 1/Ps*single(~Mask);
end
    
end