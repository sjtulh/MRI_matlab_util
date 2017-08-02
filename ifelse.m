function [res]= ifelse(flag, A, B)

if(flag)
    res = A;
else
    res = B;
end

end