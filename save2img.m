function save2img(filename, img, scale, map)

img = double(img);

if nargin < 4
    map = colormap('gray');
end

if nargin < 3
    mini = 0;
    maxi = max(abs(img(:)));
elseif length(scale) == 1
    mini = 0;
    maxi = scale;
else
    mini = scale(1);
    maxi = scale(2);
end

if ndims(squeeze(img)) > 2
    
    imwrite(img, filename);
    
else

    img(img < mini) = mini;
    img(img > maxi) = maxi;

    img = (img - mini) / (maxi - mini);

    if ~isempty(regexp(filename, '.jpe?g$'))
        imwrite(img*255, map, filename);
    else
        imwrite(map_rgb(img*63+1, map), filename);
    end

end

end

function img_rgb = map_rgb(img_gray, map)

img_rgb = reshape(map([floor(img_gray(:))],:), [size(img_gray),3]);

end