function flist = read_flist(folder_name, idx_file)

if nargin < 2
    idx_file = 0;
end

flist = cellstr(ls(folder_name));
flist = flist(3:end,:);

if idx_file > 0
    flist = flist{idx_file};
elseif idx_file < 0
    flist = flist{end+idx_file+1};
end

end