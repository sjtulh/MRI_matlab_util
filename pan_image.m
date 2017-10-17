function res = pan_image(img, flag_switch)

if nargin < 2
    flag_switch = 0;
end

height_img = size(img,1);
width_img = size(img,2);
if ~flag_switch
    n_row = size(img,3);
    n_col = size(img,4);
else
    n_row = size(img,4);
    n_col = size(img,3);
end
n_channel = size(img,5);

res = zeros([height_img*n_row, width_img*n_col, n_channel], 'single');
for i = 1:n_row
    for j = 1:n_col
        if ~flag_switch
            res((i-1)*height_img+(1:height_img), (j-1)*width_img+(1:width_img), :) = squeeze(img(:,:,i,j,:));
        else
            res((i-1)*height_img+(1:height_img), (j-1)*width_img+(1:width_img), :) = squeeze(img(:,:,j,i,:));
        end
    end
end

end