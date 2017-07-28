% Wrapped recon function of UTE 3D data
% Input:
%   pfile_name:     Pfile name with full path
%   iField_name:  Output name for iField and all parameters with full
%                   path
%   param_recon:    Options for recon
% ==== by Zhe Liu, 8/18/2014 ==== %

% Add features for original variable TE:
%   Read 1 more projection
% ==== by Zhe Liu, 9/23/2014 ==== %

% Add gridding method with fessler irt library
% ==== by Zhe Liu, 11/9/2014 ==== %

% Add extra projections for noise estimation
% ==== by Zhe Liu, 5/15/2015 ==== %

% Suitable PSD:
%   pr16_id21



function UTE_recon_simple(pfile_name, iField_name, param_recon)

if nargin < 3
    param_recon = [];
end

if nargin < 2
    iField_name = 'iField';
end

% Set default recon options
param_recon_def.n_cut_point = 0;
param_recon_def.offset_gradient = zeros(1,3,'double');
param_recon_def.flag_reg = 0;
param_recon_def.DCF_method = 'voronoi';
param_recon_def.flag_save_per_echo = 0;
param_recon_def.flag_split_echo = 0;
param_recon_def.flag_SENSE = 0;
param_recon_def.flag_noiseest = 0;
param_recon_def.TE_list = [];
param_recon_def.SE_list = [];
param_recon_def.offset_fov = [0, 0, 0];
param_recon_def.flag_combine_coil = 1;

% Load recon options
tmp = fieldnames(param_recon_def);
for i = 1:length(tmp)
    each_field = tmp{i};
    if isfield(param_recon, each_field) 
        eval([each_field, ' = getfield(param_recon, ''', each_field, ''');']);
    else
        eval([each_field, ' = getfield(param_recon_def, ''', each_field, ''');']);
        eval(['param_recon.', each_field, ' = ', each_field, ';']);
    end
end



%% Main recon process

% Load header and data
[data,dummy,dummy,dummy,dummy,hdr] = read_pfile(pfile_name);

% Load parameter from header
FOV = hdr.fov/1000;
CF = double(hdr.ps_aps_freq/10);
nslice = hdr.nslices;
necho = hdr.nechoes;
ncoil = hdr.dab_stop_rcv(1)+1;
z_scale = hdr.user2;
a_gxw = double(hdr.user4);
pw_gxwa = double(hdr.user5);
opplane = hdr.user7;
delta_t= double(hdr.user8);
rhfrsize = double(hdr.user9); 
GAM = double(hdr.user10);
opxres = hdr.user11;
num_proj = hdr.user15;
TR = double(hdr.user31) * 1e-6;

const = double(FOV/opxres*1e-4);
xres = opxres;
fov = FOV*1000;
matrix_size = double([xres,xres,xres]);
voxel_size = double(fov * 1./matrix_size).*[1,1,z_scale];
if isempty(TE_list)
    TE = double([hdr.user23, hdr.user25, hdr.user24, hdr.user26]) * 1e-6;
else
    TE = TE_list * 1e-6;
end
if isempty(SE_list)
    SE_list = 1:length(TE);
end
switch opplane
    case 1  % Axial
        B0_dir = [0,0,1]';
    case 2  % Sagittal
        B0_dir = [-1,0,0]';
    case 3  % Coronal
        B0_dir = [-1,0,0]';
    otherwise  % Oblique 
        printf('Error: Unsupported orientation\n');
        B0_dir = [0,0,1]';
end

% Re-arrange data
data = permute(data,[1 2 5 3 4]);
data = reshape(data, [rhfrsize, 512*nslice, ncoil, necho]);
data = data(n_cut_point+1:end, 1:num_proj, :,:);
data = double(data);

% Calculate trajectory with gradient offset 
%   same for each echo
kx = calc_k(rhfrsize, delta_t, pw_gxwa, a_gxw, GAM, const, offset_gradient(1));
ky = calc_k(rhfrsize, delta_t, pw_gxwa, a_gxw, GAM, const, offset_gradient(2));
kz = calc_k(rhfrsize, delta_t, pw_gxwa, a_gxw, GAM, const, offset_gradient(3));

zmod = linspace(1, -1, num_proj);
xmod = cos(sqrt(2*num_proj*3.14159).*asin(zmod)).*(1.0-zmod.*zmod).^0.5;
ymod = sin(sqrt(2*num_proj*3.14159).*asin(zmod)).*(1.0-zmod.*zmod).^0.5;    

kx_mat = bsxfun(@times, kx, xmod);
ky_mat = bsxfun(@times, ky, ymod);
kz_mat = bsxfun(@times, kz, zmod);

crds = cat(3, kx_mat, ky_mat, kz_mat);
crds = permute(crds, [3,1,2]);
crds = double(crds);

% Calculate DCF
switch DCF_method
    case 'sdc' 
        fprintf('sdc: Pipe method\n');
        
        disp('sdc: Pipe method')
        effMtx = double(opxres);
        numIter = 25;
        osf     = 2.1;
        verbose = 1;

        DCF = sdc3_MAT(crds, numIter, effMtx, verbose, osf);
        DCF = double(DCF);   
        
    case 'voronoi'
        fprintf('sdc: voronoi method\n');
       
        kmax = (kx_mat(rhfrsize,1).^2+ky_mat(rhfrsize,1).^2+kz_mat(rhfrsize,1).^2).^0.5;
        dk = 2*kmax/opxres;
        ktmp = (kx_mat(:,1).^2+ky_mat(:,1).^2+kz_mat(:,1).^2).^0.5;
        d_tmp = [0,0,0;diff(kx_mat(:,1)),diff(ky_mat(:,1)),diff(kz_mat(:,1))];
        dktmp = sum(d_tmp.^2,2).^0.5;
    
        kweight = (ktmp/kmax).^2.*dktmp/dk;
        kweight = kweight; 
        DCF = repmat(kweight,1,num_proj);
        DCF = double(DCF); 
        
        % Normalization for different number of projection        
        DCF = DCF / num_proj;           
end

% Gridding
%   using Gmri as gridding kernel
nufft_args = {matrix_size, 5*ones(size(matrix_size)), 1.5*matrix_size, matrix_size/2 + offset_fov, 'table', 2^10, 'minmax:kb'};
data = reshape(data, [rhfrsize*num_proj, ncoil, necho]);
crds = permute(reshape(crds, [3, rhfrsize*num_proj]),[2,1]);
DCF = DCF(:);

mask = true(matrix_size);
basis = {'rect'};
kspace = repmat((1./(fov./matrix_size)), [rhfrsize*num_proj, 1]).*crds;

Gm = Gmri(kspace, mask, 'fov', fov, 'basis', basis, 'nufft', nufft_args);

iField = zeros([matrix_size, 1, ncoil], 'single');

for i_echo = 1:necho
    for i_coil = 1:ncoil       
        fprintf('i_echo = %d, i_coil = %d\n', i_echo, i_coil);
        
        iField(:,:,:,1,i_coil) = reshape(Gm'*(DCF.*data(:,i_coil,i_echo)), matrix_size);
    end
    save([iField_name, '_grid_echo', num2str(i_echo), '.mat'], 'iField', 'offset_fov', '-v7.3')   
end

% Combine coil
if flag_combine_coil

    N_echo = length(SE_list);
    
    iField_sos = zeros([matrix_size, N_echo], 'single');
    iField_comp = zeros([matrix_size, N_echo], 'single');

    % Load and combine data 
    for i_SE = 1:length(SE_list)

        fprintf('i_SE=%d\n', i_SE);

        load([iField_name, '_grid_echo', num2str((i_SE)), '.mat' ], 'iField');

        iField_sos(:,:,:,i_SE) = sqrt(sum(abs(iField).^2,5));

        if i_SE == 1
            echo_1 = squeeze(iField);
        end

        iField_comp(:,:,:,i_SE) = sum(squeeze(iField).*conj(echo_1), 4);
        iField_comp(:,:,:,i_SE) = sqrt(abs(iField_comp(:,:,:,i_SE))).*exp(sqrt(-1)*angle(iField_comp(:,:,:,i_SE)));
    
    end

    save([iField_name, '_sos_grid.mat'],'iField_sos','CF','TE','matrix_size','voxel_size','B0_dir','-v7.3');


    clear iField echo_1 echo_i

    save([iField_name, '_comp_grid.mat'],'iField_comp','TE','CF','B0_dir','matrix_size', 'voxel_size','TR', 'param_recon', '-v7.3')

end



end






function [r_k] = calc_k(rhfrsize, delta_t, pw_gxwa, a_gxw, GAM, const, offset, flag_ramp_down)

if nargin < 8
    flag_ramp_down = 0;
end

if nargin < 7
    offset = 0;
end


ramp_npt = double(floor(pw_gxwa/delta_t));


r_k = zeros([1,rhfrsize], 'double');

if flag_ramp_down
    % With ramp down sampling
    plateau_npt = double(rhfrsize - 2*ramp_npt);
    km = GAM*a_gxw*(rhfrsize-ramp_npt)*delta_t;
    
    if offset < 0
        idx_turn = floor(ramp_npt+1+offset); 
        r_k(1:idx_turn) = GAM*0.5*a_gxw*((0:idx_turn-1) - offset).^2*delta_t.^2/pw_gxwa;
        r_k(idx_turn+1:idx_turn+plateau_npt) = GAM*a_gxw*(( 0.5*pw_gxwa )+((idx_turn+1:idx_turn+plateau_npt)-(ramp_npt+1+offset))*delta_t);
        r_k(idx_turn+plateau_npt+1:idx_turn+plateau_npt+ramp_npt) = km - GAM*0.5*a_gxw*(rhfrsize+1+offset - (idx_turn+plateau_npt+1:idx_turn+plateau_npt+ramp_npt)).^2*delta_t.^2/pw_gxwa;
        r_k(idx_turn+plateau_npt+ramp_npt+1:rhfrsize) = km;
    elseif offset > 0
        idx_turn = floor(1+offset);
        r_k(1:idx_turn) = 0;
        r_k(idx_turn+1:idx_turn+ramp_npt) = GAM*0.5*a_gxw*((idx_turn+1:idx_turn+ramp_npt) - (1+offset)).^2*delta_t.^2/pw_gxwa;
        r_k(idx_turn+ramp_npt+1:idx_turn+ramp_npt+plateau_npt) = GAM*a_gxw*(( 0.5*pw_gxwa )+((idx_turn+ramp_npt+1:idx_turn+ramp_npt+plateau_npt)-(ramp_npt+1+offset))*delta_t);
        r_k(idx_turn+ramp_npt+plateau_npt+1:rhfrsize) = km - GAM*0.5*a_gxw*(rhfrsize+1+offset - (idx_turn+ramp_npt+plateau_npt+1:rhfrsize)).^2*delta_t.^2/pw_gxwa;
    else
        r_k(1:ramp_npt+1)=GAM*0.5*a_gxw*(0:ramp_npt).^2*delta_t.^2/pw_gxwa;
        r_k(ramp_npt+2:rhfrsize-(ramp_npt-1))=GAM*a_gxw*(( 0.5*pw_gxwa )+( (ramp_npt+1:rhfrsize-(ramp_npt-1)-1)*delta_t-pw_gxwa));
        r_k(rhfrsize-(ramp_npt-1)+1:rhfrsize)=km - r_k(ramp_npt:-1:2); 
    end
else
    % No ramp down sampling
    if offset < 0
        idx_turn = floor(ramp_npt+1+offset);
        r_k(1:idx_turn) = GAM*0.5*a_gxw*((0:idx_turn-1) - offset).^2*delta_t.^2/pw_gxwa;
        r_k(idx_turn+1:rhfrsize) = GAM*a_gxw*(( 0.5*pw_gxwa )+((idx_turn+1:rhfrsize)-(ramp_npt+1+offset))*delta_t);
    elseif offset > 0
        idx_turn_1 = floor(offset);
        idx_turn_2 = floor(ramp_npt+1+offset);
        r_k(1:idx_turn_1+1) = 0;
        r_k(idx_turn_1+2:idx_turn_2) = GAM*0.5*a_gxw*((idx_turn_1+1:idx_turn_2-1) - offset).^2*delta_t.^2/pw_gxwa;
        r_k(idx_turn_2+1:rhfrsize) = GAM*a_gxw*(( 0.5*pw_gxwa )+((idx_turn_2+1:rhfrsize)-(ramp_npt+1+offset))*delta_t);
    else
        r_k(1:ramp_npt+1)=GAM*0.5*a_gxw*(0:ramp_npt).^2*delta_t.^2/pw_gxwa;
        r_k(ramp_npt+2:rhfrsize)=GAM*a_gxw*(( 0.5*pw_gxwa )+( (ramp_npt+1:rhfrsize-1)*delta_t-pw_gxwa));
    end
end

r_k = const* r_k.';

end




