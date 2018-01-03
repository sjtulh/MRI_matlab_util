function save2avi(data, range, filename)
% save file to avi

[na,nb,nc] = size(data);
id = double(data);

if isempty(range)
    range = [0, max(data(:))];
end

mi = double(range(1));
ma = double(range(2));
id = (id - mi) / (ma-mi);
id = reshape(id,[na,nb,1,nc]);
id = repmat(id,[1,1,3,1]);

idmov = immovie(id);

try
    v = VideoWriter(filename, 'Uncompressed AVI');
    v.FrameRate = 15;
    open(v)
    writeVideo(v, idmov);
    close(v)
    
    
catch
%     VideoWriter(idmov, filename);
end

end
