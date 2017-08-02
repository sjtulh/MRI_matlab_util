function res = unwrapPhase_crop(iMag, iFreq_raw, matrix_size, idxs_crop)

if nargin < 4
    idxs_crop = [1,size(iFreq_raw,1),1,size(iFreq_raw,2),1,size(iFreq_raw,3)];
end

res = embed_crop(unwrapPhase(trim3(iMag, idxs_crop), trim3(iFreq_raw, idxs_crop), [idxs_crop(2)-idxs_crop(1)+1, idxs_crop(4)-idxs_crop(3)+1, idxs_crop(6)-idxs_crop(5)+1]), matrix_size, idxs_crop);

end