function res = trim3(img, input_arg)

if size(input_arg,1) == 1 
    idxs_crop = input_arg;
    if size(input_arg,2) == 3
        idxs_crop = [1, input_arg(1), 1, input_arg(2), 1, input_arg(3)];
    end
elseif size(input_arg,2) == 1
    idxs_crop = input_arg;
    idxs_crop = [1, idxs_crop(2)-idxs_crop(1)+1, 1, idxs_crop(4)-idxs_crop(3)+1, 1, idxs_crop(6)-idxs_crop(5)+1];
else
    Mask = input_arg;

    idxs_crop = get_idxs_crop(Mask);
end

res = img(idxs_crop(1):idxs_crop(2), idxs_crop(3):idxs_crop(4), idxs_crop(5):idxs_crop(6), :);

end