function [img_fill] = imfill3(img, dim)

if nargin < 2
    dim = 7;    % All 3 dimensions
end

if bitand(dim, 4)
    img_tmp = img;
    for i=1:size(img, 1)
        i
        img_tmp(i,:,:) = permute(imfill(permute(img(i,:,:), [3,2,1]), 'holes'), [3,2,1]);
    end
    img = img_tmp;
end

if bitand(dim, 2)
    for i=1:size(img, 2)
        i
        img_tmp(:,i,:) = permute(imfill(permute(img(:,i,:), [1,3,2]), 'holes'), [1,3,2]);
    end
    img = img_tmp;
end

if bitand(dim, 1)
    for i=1:size(img, 3)
        i
        img_tmp(:,:,i) = permute(imfill(permute(img(:,:,i), [1,2,3]), 'holes'), [1,2,3]);
    end
    img = img_tmp;
end

img_fill = img;

end