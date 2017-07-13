function res = unwrapping_gc_crop(iFreq_raw, iMag, Mask, voxel_size, idxs_crop);

if nargin < 5
    idxs_crop = [1,size(iFreq_raw,1),1,size(iFreq_raw,2),1,size(iFreq_raw,3)];
end

matrix_size = size(iFreq_raw);

iFreq = unwrapping_gc(trim3(iFreq_raw, idxs_crop), trim3(iMag, idxs_crop), voxel_size);
Mask = trim3(Mask, idxs_crop);

tmp = iFreq;
idxs_try = [-5:5];
clear tmps
for i = 1:length(idxs_try)
    tmps(i) = abs(sum(col((tmp - 2*pi*idxs_try(i)).*Mask)));
end
[~,idx_tmp] = min(tmps);
idx_wrap = idxs_try(idx_tmp)
iFreq = tmp - 2*pi*idx_wrap;

res = embed_crop(iFreq, matrix_size, idxs_crop);

end