function [sens,FLAG,RELRES,ITER,RESVEC,LSVEC] = calc_SENSEMAP(x, nIter, fermi_map, imgref)
% Calculate 3D Sensitivity Map for SENSE reconstrution
% Note:
%   Implementation of image-domain SPIRiT reconstruction from arbitrary
%   k-space. The function is based on Jeff Fessler's nufft code and LSQR 
% Ref:
%   Uecker, M., et al. (2014). MRM 71(3): 990-1001.

N = numel(x);
imSize = [size(x), 1];
dataSize = imSize;

b = x(:);

[fsens,FLAG,RELRES,ITER,RESVEC,LSVEC] = lsqr(@(x,tflag)afun(x,imgref,dataSize,imSize,fermi_map,tflag), b, 1e-4, nIter, speye(N,N), speye(N,N), zeros(N,1));
fsens = reshape(fsens,imSize);
% sens = ifft3b(fsens.*repmat(sqrt(fermi_map),[1,1,1,imSize(4)]));
sens = ifft3b(fsens.*repmat(fermi_map,[1,1,1,imSize(4)]));

end


function [y, tflag] = afun(x,imgref,dataSize,imSize,fermi_map,tflag)
if strcmp(tflag,'transp')
    x = reshape(x(1:prod(dataSize)),imSize);
    y = repmat(imgref,[1,1,1,imSize(4)]).*x;
    y = repmat(sqrt(fermi_map),[1,1,1,imSize(4)]).*fft3b(y);
    y = y(:);
else
    x = reshape(x,imSize);
    x1 = ifft3b(x.*repmat(sqrt(fermi_map),[1,1,1,imSize(4)]));
    y = repmat(imgref,[1,1,1,imSize(4)]).*x1;
    y = y(:);
end
end


function fimg = fft3b(img)
% Forward 3D fft with balanced weight 1/sqrt(N)
%   (as supposed to Matlab's 1 weight on fft and 1/N on ifft)
    N = size(img,1)*size(img,2)*size(img,3);
    fimg = fftshift(fft(ifftshift(img,1),[],1),1);
    fimg = fftshift(fft(ifftshift(fimg,2),[],2),2);
    fimg = 1/sqrt(N)*fftshift(fft(ifftshift(fimg,3),[],3),3);
end


function img = ifft3b(fimg)
% Inverse 3D fft with balanced weight 1/sqrt(N)
%   (as supposed to Matlab's 1 weight on fft and 1/N on ifft)
    N = size(fimg,1)*size(fimg,2)*size(fimg,3);
    img = fftshift(ifft(ifftshift(fimg,1),[],1),1);
    img = fftshift(ifft(ifftshift(img,2),[],2),2);
    img = sqrt(N)*fftshift(ifft(ifftshift(img,3),[],3),3);
end
