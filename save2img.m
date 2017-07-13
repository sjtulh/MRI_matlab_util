function save2img(filename, img, scale, map)

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

    imwrite(img*256, map, filename);

end

end