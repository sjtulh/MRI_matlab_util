function img = load4nii(filename, flag_noflip)

if nargin < 2
    flag_noflip = 0;
end

nii = load_nii(filename);

if flag_noflip
    img = nii.img;  
else
    img = flipdim(flipdim(flipdim(nii.img,3),2),1);
end


end
