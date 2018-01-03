function res = kpad(img, fac)
matrix_size = size(img);
matrix_size = matrix_size(1:3);
matrix_size_new = matrix_size .* fac;

kspace = fftshift(fftshift(fftshift(fft(fft(fft(ifftshift(ifftshift(ifftshift(img, 1), 2), 3), [], 1), [], 2), [], 3), 1), 2), 3);

kspace = embed_crop(kspace, [matrix_size_new, size(img,4)], [1, matrix_size(1), 1, matrix_size(2), 1, matrix_size(3)]);
kspace = circshift(kspace, -[floor(matrix_size/2), 0]);

res = fftshift(fftshift(fftshift(ifft(ifft(ifft(kspace, [], 1), [], 2), [], 3), 1), 2), 3);

end