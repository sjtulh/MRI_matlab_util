function res = ktrim3c(img, dims_new)

kspace = fftshift(fftshift(fftshift(fft(fft(fft(ifftshift(ifftshift(ifftshift(img, 1), 2), 3), [], 1), [], 2), [], 3), 1), 2), 3);

kspace = trim3c(kspace, dims_new);

res = fftshift(fftshift(fftshift(ifft(ifft(ifft(ifftshift(ifftshift(ifftshift(kspace, 1), 2), 3), [], 1), [], 2), [], 3), 1), 2), 3);

end