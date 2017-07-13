function res = pick_major(Mask, n_major)

if nargin < 2
    n_major = 1;
end

CC = bwconncomp(Mask,6);
numPixels = cellfun(@numel,CC.PixelIdxList);
[numPixels_sorted,idxs] = sort(numPixels,2,'descend');
res = zeros(size(Mask));
for i = 1:n_major
    numPixels_sorted(i)
    idx = idxs(i);
    res(CC.PixelIdxList{idx}) = 1;
end

end