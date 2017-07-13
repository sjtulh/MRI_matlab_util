function res = img_window(img, window)

mini = window(1);
maxi = window(2);

img(img < mini) = mini;
img(img > maxi) = maxi;

res = (img - mini) / (maxi - mini);

end