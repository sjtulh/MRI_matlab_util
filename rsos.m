function res = rsos(img)

res = sqrt(sum(abs(img).^2, 4));

end