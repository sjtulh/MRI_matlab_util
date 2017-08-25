function res = dilate3_SMV(img, voxel_size, radius)

res = SMV(single(img), size(img), voxel_size, radius) > 0.001;

end