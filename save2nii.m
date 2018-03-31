function save2nii(img, voxel_size, filename)

nii = make_nii(flipdim(flipdim(single(img),1),2),voxel_size,[0,0,0],16);     % 16 for float32
% nii = make_nii(single(img),voxel_size,[0,0,0],16);     % 16 for float32
save_nii(nii, filename);

end
