function res = trim3(img, input_arg)

if size(input_arg,1) == 1 
    idxs_crop = input_arg;
elseif size(input_arg,2) == 1
    idxs_crop = input_arg;
    idxs_crop = [1, idxs_crop(2)-idxs_crop(1)+1, 1, idxs_crop(4)-idxs_crop(3)+1, 1, idxs_crop(6)-idxs_crop(5)+1];
else
    Mask = input_arg;

    matrix_size0 = size(img);
    d1 = max(max(Mask,[],2),[],3);
    d1first = find(d1,1,'first');
    d1last = find(d1,1,'last');

    d2 = max(max(Mask,[],1),[],3);
    d2first = find(d2,1,'first');
    d2last = find(d2,1,'last');

    d3 = max(max(Mask,[],1),[],2);
    d3first = find(d3,1,'first');
    d3last = find(d3,1,'last');

    idxs_crop = [d1first, d1last, d2first, d2last, d3first, d3last]
end

res = img(idxs_crop(1):idxs_crop(2), idxs_crop(3):idxs_crop(4), idxs_crop(5):idxs_crop(6), :);

end