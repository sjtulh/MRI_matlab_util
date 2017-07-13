function res = embed_crop(img, matrix_size_raw, idxs_crop, val_fill)

if nargin < 4
    val_fill = 0;
end

flag_logical = islogical(img);
if flag_logical
    img = single(img);
end

if size(idxs_crop,2) == 1
    idxs_crop = [1, idxs_crop(2)-idxs_crop(1)+1, 1, idxs_crop(4)-idxs_crop(3)+1, 1, idxs_crop(6)-idxs_crop(5)+1];
end

res = ones(matrix_size_raw, class(img))*val_fill;
res(idxs_crop(1):idxs_crop(2), idxs_crop(3):idxs_crop(4), idxs_crop(5):idxs_crop(6),:) = img(1:idxs_crop(2)-idxs_crop(1)+1, 1:idxs_crop(4)-idxs_crop(3)+1, 1:idxs_crop(6)-idxs_crop(5)+1,:);

if flag_logical
    res = logical(res);
end

end