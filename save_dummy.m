function save_dummy(fn, n_item, varargin)

vararg_save = {};
vararg_save{1} = fn;
% str_cmd = ['save ', fn];
idx = 1;
for i_item = 1:n_item
    eval([varargin{idx}, ' = varargin{', num2str(idx+1),'};']);
    vararg_save{i_item+1} = varargin{idx};
%     str_cmd = [str_cmd, ' ', varargin{idx}];
    idx = idx+2;
end
vararg_save = [vararg_save, varargin(idx:end)];
% vararg_save{:}
% str_cmd = [str_cmd, ' ', varargin{idx:end}]

save(vararg_save{:});

end