% Compute approximate distance map (Euclidean) using Matlab built-in bwdist.m
%   with anisotropic voxel size using zero-padding
%   Trick:
%       1. Use zero-padding to achieve isotropic voxel size in order to 
%           use bwdist.m
%       2. 
%           Option 2.1: Linear interpolation
%           Option 2.2: K-space cropping
%

function d_map = bwdistzp(img, voxel_size)

matrix_size = size(img);
mx = matrix_size(1);
my = matrix_size(2);
mz = matrix_size(3);

% Anisotropic voxel size [vx, vy, vz]
%   Assume: vx = vy < vz
%       Only need to zero-pad in z-direction
vx = voxel_size(1);
vy = voxel_size(2);
vz = voxel_size(3);

% fx = mx*vx;
% fy = my*vy;
% fz = mz*vz;

mz_pad = floor(matrix_size(3)*vz/vx);
matrix_size_pad = [mx, my, mz_pad];

% Zero-padding in k-space
shift_z = floor(mz/2);
fimg = fftn(ifftshift(img));
fimg_pad = circshift(padarray(circshift(fimg, [0,0,shift_z]), double([0,0,mz_pad-mz]), 0, 'post'), [0,0,-shift_z]);
img_pad = fftshift(ifftn(fimg_pad));
img_pad = img_pad > 0.1;

% Calculate dist map on padded imag
d_map_pad = bwdist(img_pad)*vx;

% % Option 2.1: Interpolation back to original voxel size
% [X_pad,Y_pad,Z_pad] = meshgrid(linspacec(0,my,my), linspacec(0,mx,mx), linspacec(0,mz_pad,mz_pad));
% [X,Y,Z] = meshgrid(linspacec(0,my,my), linspacec(0,mx,mx), linspacec(0,mz_pad,mz));
% d_map = interp3(X_pad, Y_pad, Z_pad, d_map_pad, X, Y, Z, 'linear');

% Option 2.2: K-space cropping
fd_map_pad = fftn(ifftshift(d_map_pad));
fd_map = circshift(trim3(circshift(fd_map_pad, [0,0,shift_z]), matrix_size), [0,0,-shift_z]);
d_map = fftshift(ifftn(fd_map));
d_map = real(d_map);
d_map(d_map < 0) = 0;

end