function res = correct_2pi_shift(iFreq, Mask)

tmp = mean(iFreq(Mask > 0));
res = iFreq - 2*pi*round(tmp/(2*pi));

end