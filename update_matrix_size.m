function [matrix_size, margin, padsize] = update_matrix_size(matrix_size, idxs_crop)

margin = max(floor(matrix_size/256*40), 40);
matrix_size = [ floor((idxs_crop(2) - idxs_crop(1)+margin(1))/2)*2,...
                floor((idxs_crop(4) - idxs_crop(3)+margin(2))/2)*2,...
                floor((idxs_crop(6) - idxs_crop(5)+margin(3))/2)*2];
padsize = double([matrix_size(1)-(idxs_crop(2) - idxs_crop(1)+1),matrix_size(2)-(idxs_crop(4) - idxs_crop(3)+1),matrix_size(3)-(idxs_crop(6) - idxs_crop(5)+1)]);

end