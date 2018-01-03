function [iFreq_corrected, n_2pi] = correct_2pi_shift(iFreq, Mask)

n_2pi = round(mean(iFreq(Mask > 0))/(2*pi));
iFreq_corrected = iFreq - 2*pi*n_2pi;

end