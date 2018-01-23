% Wrapped recon function of GE data
% Input:
%   pfile_name:     Pfile name with full path
%   iField_name:  Output name for iField and all parameters with full
%                   path
%   param_recon:    Options for recon
% ==== by Zhe Liu, 1/14/2018 ==== %
%
% Suitable PSD:
%   me



function ME_recon(pfile_name, iField_name, param_recon)

if nargin < 3
    param_recon = [];
end

if nargin < 2
    iField_name = 'iField';
end

% Set default recon options
param_recon_def = ME_default_param();
  
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

% -------------------- (Work for DV26) --------------------- %
[data, hdr] = read_pfile(pfile_name);
data = permute(data, [1,2,4,3,5]);
% Load parameter from header
FOV = double(hdr.rdb_hdr_image.dfov/1000);
CF = double(hdr.rdb_hdr_ps.aps_freq/10);
nslice = double(hdr.rdb_hdr_rec.rdb_hdr_nslices);
necho = double(hdr.rdb_hdr_rec.rdb_hdr_nechoes);
ncoil = double(hdr.rdb_hdr_rec.rdb_hdr_dab.stop_rcv(1))-...
    double(hdr.rdb_hdr_rec.rdb_hdr_dab.start_rcv(1))+1;
% opplane = double(hdr.rdb_hdr_rec.rdb_hdr_user7);
if scan_orient == 'A'
    opplane = 1;
elseif scan_orient == 'S'
    opplane = 2;
elseif scan_orient == 'C'
    opplane = 3;
end
if isempty(TE_list)
    TE = double([hdr.rdb_hdr_rec.rdb_hdr_user23, hdr.rdb_hdr_rec.rdb_hdr_user25, hdr.rdb_hdr_rec.rdb_hdr_user24, hdr.rdb_hdr_rec.rdb_hdr_user26]) * 1e-6;
else
    TE = TE_list * 1e-6;
end

fov = FOV*1000;
matrix_size = double([size(data,1), size(data,2), size(data,3)]);
% voxel_size = double(fov * 1./matrix_size);

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
% --------- Extra pre-processing ------- %
% 1. Select coil (channel)
if ~isempty(coil_select)
    data = data(:, :, :, :, coil_select);
    ncoil = length(coil_select);
end    
data = double(data);

img = fftshift(fftshift(ifft(ifft(ifft(ifftshift(ifftshift(ifftshift(data,1),2),3),[],1),[],2),[],3),1),2);


%% Combine echo
if flag_combine_echo
    
    N_echo = length(SE_list);

    iField_sos = zeros([matrix_size, N_echo], 'single');
    iField_comp = zeros([matrix_size, N_echo], 'single');  
        
    for i_SE = 1:length(SE_list)

        fprintf('i_SE=%d\n', i_SE);

        if ~flag_SENSE
            % Gridding, need for coil combination
            iField = squeeze(img(:,:,:,i_SE,:));
            
            iField_sos(:,:,:,i_SE) = sqrt(sum(abs(iField).^2, 4));

            if i_SE == 1
                echo_1 = squeeze(iField);
            end

            iField_comp(:,:,:,i_SE) = sum(squeeze(iField).*conj(echo_1), 4);
            iField_comp(:,:,:,i_SE) = sqrt(abs(iField_comp(:,:,:,i_SE))).*exp(sqrt(-1)*angle(iField_comp(:,:,:,i_SE)));
        else
            % SENSE, no need for coil combination
            load([iField_name, '_SENSE_echo', num2str((i_SE)), '.mat' ], 'iField');

            iField_comp(:,:,:,i_SE) = iField;
            iField_sos(:,:,:,i_SE) = abs(iField);
        end

    end

    save([iField_name, '_sos.mat'],'iField_sos','CF','TE','matrix_size','voxel_size','B0_dir','-v7.3');


    clear iField echo_1 echo_i

    save([iField_name, '_comp.mat'],'iField_comp','TE','CF','B0_dir','matrix_size', 'voxel_size', 'param_recon', '-v7.3')

end



end




