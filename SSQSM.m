function [chi_ss, mask_sharp] = wrapper_SSQSM(magn, phs, mask_bet, voxel_size, TE, B0, B0_dir, gyro, param)

%% =========== from Bilgic's script_3D_EPI_SSQSM =======================
if nargin < 9
    param = [];
end

N = size(mask_bet);
voxel_size = col(voxel_size)';

lambda = 2.9e-2;
if isfield(param, 'lambda')
    lambda = param.lambda;
end

tol_pcg = 1e-2;
if isfield(param, 'tol_pcg')
    tol_pcg = param.tol_pcg;
end

maxiter_pcg = 100;
if isfield(param, 'maxiter_pcg')
    maxiter_pcg = param.maxiter_pcg;
end

diameter_in = 11;
if isfield(param, 'diameter_in')
    diameter_in = param.diameter_in;
end

diameter_out = 3;
if isfield(param, 'diameter_out')
    diameter_out = param.diameter_out;
end

step_diameter = 2;
if isfield(param, 'step_diameter')
    step_diameter = param.step_diameter;
end



%% Laplacian unwrapping kernel

ksize = [3, 3, 3];               
khsize = (ksize-1)/2;

kernel = [];

kernel(:,:,1) = [0 0 0; 0 1 0; 0 0 0];
kernel(:,:,2) = [0 1 0; 1 -6 1; 0 1 0];
kernel(:,:,3) = [0 0 0; 0 1 0; 0 0 0];

Kernel = zeros(N);
Kernel( 1+N(1)/2 - khsize(1) : 1+N(1)/2 + khsize(1), 1+N(2)/2 - khsize(2) : 1+N(2)/2 + khsize(2), 1+N(3)/2 - khsize(3) : 1+N(3)/2 + khsize(3) ) = -kernel;

del_op = fftn(fftshift(Kernel));
del_inv = zeros(size(del_op));
del_inv( del_op~=0 ) = 1 ./ del_op( del_op~=0 );

exp_phase = exp(1i.*phs);
Phase_unwrap = del_inv .* fftn(imag(conj(exp_phase) .* ifftn(del_op .* fftn(exp_phase))));

% For highly unisotropic voxelsize
%   input iFreq rather than iFreq_raw
% ==== By Zhe Liu, 3/11/2016 ==== %
Phase_unwrap = fftn(phs);
phs_bk = phs;

%% V-Sharp filtering for background removal
stol = .2;                  % truncation threshold
Kernel_Sizes = diameter_in:-step_diameter:diameter_out;
% Kernel_Sizes = 3;           % Changed from Bilgic's script due to small matrix size

DiffMask = zeros([N, length(Kernel_Sizes)]);
Mask_Sharp = zeros([N, length(Kernel_Sizes)]);
Del_Sharp = zeros([N, length(Kernel_Sizes)]);

phs = 0;

% for k = 1:length(Kernel_Sizes)
%     Kernel_Size = Kernel_Sizes(k);
% 
%     ksize = [Kernel_Size, Kernel_Size, Kernel_Size];                % Sharp kernel size
% 
%     khsize = (ksize-1)/2;
%     [a,b,c] = meshgrid(-khsize(2):khsize(2), -khsize(1):khsize(1), -khsize(3):khsize(3));
% 
%     kernel = (a.^2 / khsize(1)^2 + b.^2 / khsize(2)^2 + c.^2 / khsize(3)^2 ) <= 1;
%     kernel = -kernel / sum(kernel(:));
%     kernel(khsize(1)+1,khsize(2)+1,khsize(3)+1) = 1 + kernel(khsize(1)+1,khsize(2)+1,khsize(3)+1);
% 
%     Kernel = zeros(N);
%     Kernel( 1+N(1)/2 - khsize(1) : 1+N(1)/2 + khsize(1), 1+N(2)/2 - khsize(2) : 1+N(2)/2 + khsize(2), 1+N(3)/2 - khsize(3) : 1+N(3)/2 + khsize(3) ) = kernel;
% 
%     del_sharp = fftn(fftshift(Kernel));    
% 
%     % erode mask to remove convolution artifacts
%     erode_size = ksize + 1;
% 
%     mask_bet = circshift(mask_bet, [1,1,1]);
%     msk_sharp = imerode(mask_bet, strel('line', erode_size(1), 0));
%     msk_sharp = imerode(msk_sharp, strel('line', erode_size(2), 90));
%     msk_sharp = permute(msk_sharp, [1,3,2]);
%     msk_sharp = imerode(msk_sharp, strel('line', erode_size(3), 0));
%     msk_sharp = permute(msk_sharp, [1,3,2]);
%     mask_bet = circshift(mask_bet, [-1,-1,-1]);
%     msk_sharp = circshift(msk_sharp, [-1,-1,-1]);    
% 
%     Mask_Sharp(:,:,:,k) = msk_sharp; 
%     Del_Sharp(:,:,:,k) = del_sharp; 
%     
%     if k == 1
%         DiffMask(:,:,:,1) = Mask_Sharp(:,:,:,1);
%     else
%         DiffMask(:,:,:,k) = Mask_Sharp(:,:,:,k) - Mask_Sharp(:,:,:,k-1);
%     end
%     
%     
%     % Use Tian's SMV kernel instead
%     del_sharp = SMV_kernel(N, voxel_size, Kernel_Size);
%     Del_Sharp(:,:,:,k) = del_sharp; 
%     
%     
%     
%     phs = phs + DiffMask(:,:,:,k) .* ifftn(Del_Sharp(:,:,:,k) .* Phase_unwrap) / (TE * B0 * gyro);
% end

% A cleaner version
[Del_Sharp, dummy, DiffMask, Mask_Sharp] = prepare_VSHARP(Kernel_Sizes/2, mask_bet, N, voxel_size);
phs = 0;
for k = 1:size(DiffMask,4)
    phs = phs + DiffMask(:,:,:,k) .* real(ifftn(Del_Sharp(:,:,:,k) .* Phase_unwrap)) / (TE * B0 * gyro);
end

del_sharp = Del_Sharp(:,:,:,1);         % first kernel
mask_sharp = Mask_Sharp(:,:,:,end);     % largest mask

inv_sharp = zeros(size(del_sharp));
inv_sharp( abs(del_sharp)>stol ) = 1 ./ del_sharp( abs(del_sharp)>stol );

phase_sharp = mask_sharp .* ifftn(inv_sharp .* fftn(phs));





%% gradient masks from magnitude image using k-space gradients

[k2,k1,k3] = meshgrid(0:N(2)-1, 0:N(1)-1, 0:N(3)-1);
fdx = -1 + exp(2*pi*1i*k1/N(1));
fdy = -1 + exp(2*pi*1i*k2/N(2));
fdz = -1 + exp(2*pi*1i*k3/N(3));

cfdx = conj(fdx);       cfdy = conj(fdy);       cfdz = conj(fdz);

magni = magn .* mask_sharp;

Magn = fftn(magni / max(magni(:)));
magn_grad = cat(4, ifftn(Magn.*fdx), ifftn(Magn.*fdy), ifftn(Magn.*fdz));

magn_weight = zeros(size(magn_grad));

for s = 1:size(magn_grad,4)
    magn_use = abs(magn_grad(:,:,:,s));
    
    magn_order = sort(magn_use(mask_sharp==1), 'descend');
    magn_threshold = magn_order( round(length(magn_order) * .15) );
%     magn_threshold = magn_order( round(length(magn_order) * .2) );  % Keep same with TFI (and MEDI)
    magn_weight(:,:,:,s) = magn_use <= magn_threshold;
end


% To be fair, set magn_weight as that in wG in TFI (or MEDI)
if isfield(param, 'magn_weight')
    magn_weight = param.magn_weight;
end





%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
%% Single-step QSM with V-Sharp
%% ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

FOV = N .* voxel_size;   % Field of view(in milimeters)
% D = fftshift(kspace_kernel(FOV, N, B0_dir));        % k-space kernel
D = dipole_kernel(N, voxel_size, B0_dir);
E2 = abs(fdx).^2 + abs(fdy).^2 + abs(fdz).^2;

Rhs_pcg = 0;
for k = 1:size(Del_Sharp,4)
    Rhs_pcg = Rhs_pcg + conj(Del_Sharp(:,:,:,k)) .* fftn(DiffMask(:,:,:,k) .* phs);
end
Rhs_pcg = conj(D) .* Rhs_pcg;       % right hand side


% lambda_ss = 2.9e-2;                 % regularization parameter
lambda_ss = lambda;                  

B_inv = 1 ./ (eps + abs(Del_Sharp(:,:,:,1) .* D).^2 + lambda_ss*E2);        % preconditioner
precond_inv = @(x, B_inv) B_inv(:).*x;

x0 = B_inv .* Rhs_pcg;              % initial guess
x0(:) = 0;


% tic
%     [F_chi, flag, pcg_res, pcg_iter] = pcg(@(x)SSQSM_vsharp(x, D, Del_Sharp, conj(Del_Sharp), DiffMask, lambda_ss, fdx, fdy, fdz, cfdx, cfdy, cfdz, magn_weight), ...
%         Rhs_pcg(:), 1e-2, 30, @(x) precond_inv(x, B_inv), [], x0(:));
    [F_chi, flag, pcg_res, pcg_iter] = pcg(@(x)SSQSM_vsharp(x, D, Del_Sharp, conj(Del_Sharp), DiffMask, lambda_ss, fdx, fdy, fdz, cfdx, cfdy, cfdz, magn_weight), ...
        Rhs_pcg(:), tol_pcg, maxiter_pcg, @(x) precond_inv(x, B_inv), [], x0(:));       % Keep same with TFI (and MEDI): 1 Newton iteration
% toc

disp(['PCG iter: ', num2str(pcg_iter), '   PCG residual: ', num2str(pcg_res)])

chi_ss = real(ifftn(reshape(F_chi, N))) .* mask_sharp;

