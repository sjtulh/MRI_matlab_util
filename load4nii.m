function [img, voxel_size] = load4nii(filename, flag_noflip)

if nargin < 2
    flag_noflip = 0;
end

nii = load_nii(filename);
voxel_size = nii.hdr.dime.pixdim;
voxel_size = voxel_size(2:4);

if flag_noflip
    img = nii.img;  
else
    img = flipdim(flipdim(nii.img,2),1);
end


end
