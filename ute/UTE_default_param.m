function param_recon_def = UTE_default_param()

% Set default recon options
param_recon_def.n_cut_point = 0;
param_recon_def.offset_gradient = zeros(1,3,'double');
param_recon_def.flag_reg = 0;
param_recon_def.DCF_method = 'voronoi';
param_recon_def.flag_save_per_echo = 1;
param_recon_def.flag_split_echo = 0;
param_recon_def.flag_SENSE = 0;
param_recon_def.flag_noiseest = 0;
param_recon_def.TE_list = [];
param_recon_def.SE_list = [];
param_recon_def.offset_fov = [0, 0, 0];
param_recon_def.flag_combine_coil = 1;
param_recon_def.flag_combine_echo = 1;
param_recon_def.fac_under_manual = 1;

% Gridding option
param_recon_def.n_neighbor_grid = 5;
param_recon_def.fac_fft_size_grid = 1.5;

% SENSE option
param_recon_def.flag_save_inter = 1;
param_recon_def.flag_use_cache = 0;
% Filter (fermi)
param_recon_def.sens_cutoff = 6;
param_recon_def.sens_transwidth = 3;
% Reg option
param_recon_def.lambda = 5e-2;
param_recon_def.cg_max_iter = 10;
param_recon_def.cg_tol = 1e-4;
param_recon_def.cg_tol_noise = 1e-10;
param_recon_def.grad = @fgrad;
param_recon_def.div = @bdiv;
param_recon_def.e_grad = 1e-6;
param_recon_def.flag_skip_LS = 0;

% Computing option
param_recon_def.flag_parallel = 0;
param_recon_def.n_parallel = 3;

end