function [Mask, thresh] = masking_auto(img_ref, flag_method)

debug = 1;

n_bin = 1000;

if nargin < 2
    flag_method = 1;
end


samples = img_ref(img_ref~=0);
N = numel(samples);
[freqs, bins] = hist(samples, n_bin);

switch flag_method
    case 1
        n_state = 3;
        [map_state, prob_state_raw, alphas_state, mus_state, sigmas_state] = EM_gaussian_mix(double(img_ref), n_state);
        ratios_change = diff(log(mus_state));
        [~, idx_switch] = max(ratios_change);
        prob_state = cat(4, sum(prob_state_raw(:,:,:,1:idx_switch),4), sum(prob_state_raw(:,:,:,idx_switch+1:end),4));
        Mask = prob_state(:,:,:,2) > 0.9999;   % Almost sure it's foreground
    
        thresh = min(img_ref(Mask));
        
    case 2
        n_state = 2;
        [map_state, prob_state_raw, alphas_state, mus_state, sigmas_state] = EM_gaussian_mix(double(img_ref), n_state);
        mask_extract = bins > mus_state(1) & bins < mus_state(2);
        freqs_extract = freqs(mask_extract);
        bins_extract = bins(mask_extract);
%         [dummy, idx_select] = min(flipdim(freqs_extract,2));
%         thresh = bins_extract(end-idx_select+1);
        idx_select = find(freqs_extract == min(freqs_extract), 1, 'first');
        thresh = bins_extract(idx_select);
        
        Mask = img_ref > thresh; 
        
    case 3  % Mixture of method 1 and 2
        n_state = 3;
        [map_state, prob_state_raw, alphas_state, mus_state, sigmas_state] = EM_gaussian_mix(double(img_ref), n_state);
        ratios_change = diff(log(mus_state));
        [~, idx_switch] = max(ratios_change);
        mask_extract = bins > mus_state(idx_switch) & bins < mus_state(idx_switch+1);
        freqs_extract = freqs(mask_extract);
        bins_extract = bins(mask_extract);
        idx_select = find(freqs_extract == min(freqs_extract), 1, 'first');
        thresh = bins_extract(idx_select);        
        
        Mask = img_ref > thresh;      
end

if debug
    figure
    [freqs, bins] = hist(img_ref(img_ref~=0), 1000);
    hist(img_ref(img_ref~=0), 1000);
    hold on
    ys = linspace(0,max(freqs(:)),1000);
    plot(ones(size(ys))*thresh, ys, 'r-');
end

end