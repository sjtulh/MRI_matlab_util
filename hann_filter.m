function  H = hann_filter(matrix_size, voxel_size, fc)   

    if numel(matrix_size)==3
        voxel_size = voxel_size/voxel_size(2);
        sy=matrix_size(2)/matrix_size(2)/voxel_size(2);
        sx=matrix_size(2)/matrix_size(1)/voxel_size(1);
        sz=matrix_size(2)/matrix_size(3)/voxel_size(3);

        H = zeros(matrix_size);

        [Y,X,Z]=meshgrid( -matrix_size(2)/2*sy:sy:sy*(matrix_size(2)/2-1), -matrix_size(1)/2*sx:sx:(matrix_size(1)/2-1)*sx, -matrix_size(3)/2*sz:sz:(matrix_size(3)/2-1)*sz);

        n = pi*sqrt(Y.^2+X.^2+Z.^2)/(fc/2);
        n(n>pi) = pi;

        W = 0.5*(1+cos(n));
        H = W;
    else
        voxel_size = voxel_size/voxel_size(2);
        sy=matrix_size(2)/matrix_size(2)/voxel_size(2);
        sx=matrix_size(2)/matrix_size(1)/voxel_size(1);

        H = zeros(matrix_size);

        [Y,X]=meshgrid( -matrix_size(2)/2*sy:sy:sy*(matrix_size(2)/2-1), -matrix_size(1)/2*sx:sx:(matrix_size(1)/2-1)*sx);

        n = pi*sqrt(Y.^2+X.^2)/(fc/2);
        n(n>pi) = pi;

        W = 0.5*(1+cos(n));
        H = W;
    end
    