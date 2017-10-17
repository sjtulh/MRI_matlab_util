function x = TKD(RDF, threshold, matrix_size, voxel_size, B0_dir)

D = dipole_kernel(matrix_size, voxel_size, B0_dir);
D1 = D;
S = sign(D1);
S(S==0) = 1;
D1(abs(D)<threshold) = S(abs(D)<threshold)*threshold;

X = fftn(RDF)./D1;

% x = ifftn(X)/(2*pi*delta_TE*CF*1e-6);
x = ifftn(X);

end