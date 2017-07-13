% Note:
%   Length should be EVEN number along each dimension

function [idxs_crop] = get_idxs_crop(Mask, margin)

if nargin < 2
    margin = 0;
end

matrix_size = size(Mask);

d1 = max(max(Mask,[],2),[],3);
d1first = max(find(d1,1,'first')-margin, 1);
d1last = min(find(d1,1,'last')+margin, matrix_size(1));

d2 = max(max(Mask,[],1),[],3);
d2first = max(find(d2,1,'first')-margin, 1);
d2last = min(find(d2,1,'last')+margin, matrix_size(2));
        
d3 = max(max(Mask,[],1),[],2);
d3first = max(find(d3,1,'first')-margin, 1);
d3last = min(find(d3,1,'last')+margin, matrix_size(3));

% Check if odd dimension
if mod(d1last - d1first + 1, 2) == 1
    if d1first > 1
        d1first = d1first - 1;
    elseif d1last < matrix_size(1)
        d1last = d1last + 1;
    else
        error('Dimension 1 must have even length.');
    end
end
if mod(d2last - d2first + 1, 2) == 1
    if d2first > 1
        d2first = d2first - 1;
    elseif d2last < matrix_size(2)
        d2last = d2last + 1;
    else
        error('Dimension 2 must have even length.');
    end
end
if mod(d3last - d3first + 1, 2) == 1
    if d3first > 1
        d3first = d3first - 1;
    elseif d3last < matrix_size(3)
        d3last = d3last + 1;
    else
        error('Dimension 3 must have even length.');
    end
end

idxs_crop = [d1first, d1last, d2first, d2last, d3first, d3last];

end