function res = erode3_SMV(img, voxel_size, radius)

res = SMV(single(img), size(img), voxel_size, radius) > 0.999;

end