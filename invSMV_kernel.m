function invL = invSMV_kernel(matrix_size, voxel_size, radius)

L = SMV_kernel(matrix_size, voxel_size, radius);
cutoff=4; index = 1;
H = fftshift(hann_filter(matrix_size, voxel_size, cutoff*index));
invL = (1./L).*(1-H);
while ((max(abs(invL(:))))>1000)
    index = index+1;
    H = fftshift(hann_filter(matrix_size, voxel_size, cutoff*index));
    invL = (1./L).*(1-H);
end
invL(isnan(invL)) = 0;
invL(isinf(invL)) = 0;
    
end