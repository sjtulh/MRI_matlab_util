function res = unwrapping_gc_crop(iFreq_raw,iMag,voxel_size,idxs_crop,SUBSAMPLE)

if nargin < 4
    idxs_crop = [1,size(iFreq_raw,1),1,size(iFreq_raw,2),1,size(iFreq_raw,3)];
end

if nargin < 5
    SUBSAMPLE = 1;
end

matrix_size = size(iFreq_raw);

res = embed_crop(unwrapping_gc(trim3(iFreq_raw, idxs_crop), trim3(iMag, idxs_crop), voxel_size, SUBSAMPLE), matrix_size, idxs_crop);

end