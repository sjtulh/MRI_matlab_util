function res = unwrapPhase_crop(iMag, iFreq_raw, idxs_crop);

if nargin < 3
    idxs_crop = [1,size(iFreq_raw,1),1,size(iFreq_raw,2),1,size(iFreq_raw,3)];
end

matrix_size = size(iFreq_raw);

iFreq = unwrapPhase(trim3(iMag, idxs_crop), trim3(iFreq_raw, idxs_crop), [idxs_crop(2)-idxs_crop(1)+1, idxs_crop(4)-idxs_crop(3)+1, idxs_crop(6)-idxs_crop(5)+1]);
res = embed_crop(iFreq, matrix_size, idxs_crop);

end