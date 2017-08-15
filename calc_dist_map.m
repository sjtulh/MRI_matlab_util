function res = calc_dist_map(Mask, voxel_size, type_dist)
% Erosion/Dilation using SMV

if nargin < 3
    type_dist = 1;
end

if nargin < 2
    voxel_size = [1,1,1];
end

Mask_raw = single(Mask);
Mask_now = single(Mask);
d = 0;
d_map = zeros(size(Mask));

switch type_dist
    case 1
        % Manhattan (block) distance
        disp('Generate distance map using Manhattan Dist.')
        
        % assume voxel_size is [1,1,1]
        
        % Dilation by 1 voxel (radius)
        d_map(Mask_now>0) = d;
        n_unset = sum(Mask_now(:)==0);
        
        while n_unset > 0
            d = d+1
            Mask_next = single(SMV(Mask_now, size(Mask), voxel_size, 1.01) > 0.001);
            d_map(Mask_next & ~Mask_now) = d;
            Mask_now = Mask_next;
            n_unset = sum(Mask_now(:)==0);
        end
end



res = d_map;

end