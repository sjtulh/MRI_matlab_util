function y = linspacec(d1, d2, n)
% Similar to linspace but y is on the CENTER of each of n bins

y = linspace(d1, d2, 2*n+1);
y = y(2:2:end);

end