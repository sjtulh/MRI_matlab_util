function res = make_clean(data, val_clean)

if nargin < 2
    val_clean = 0;
end

res = data;
res(isnan(res) | isinf(res)) = val_clean;

end