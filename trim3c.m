function res = trim3c(img, dims_new)

matrix_size = size(img);
matrix_size = matrix_size(1:3);

idxs_crop = [(matrix_size(1)/2-dims_new(1)/2)+[1,dims_new(1)], (matrix_size(2)/2-dims_new(2)/2)+[1,dims_new(2)], (matrix_size(3)/2-dims_new(3)/2)+[1,dims_new(3)]];

res = trim3(img, idxs_crop);

end