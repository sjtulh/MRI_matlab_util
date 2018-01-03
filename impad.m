function res = impad(img, fac)
matrix_size = size(img);
matrix_size = matrix_size(1:3);
matrix_size_new = matrix_size .* fac;
[X, Y, Z] = meshgrid(linspace(0, 1, matrix_size(1)), linspace(0, 1, matrix_size(2)), linspace(0, 1, matrix_size(3)));
[X_new, Y_new, Z_new] = meshgrid(linspace(0, 1, matrix_size_new(1)), linspace(0, 1, matrix_size_new(2)), linspace(0, 1, matrix_size_new(3)));

res = interp3(X, Y, Z, img, X_new, Y_new, Z_new, 'nearest');

end