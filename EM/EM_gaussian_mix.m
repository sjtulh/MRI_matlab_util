function [states, prob_states, alphas_state, mus_state, sigmas_state] = EM_gaussian_mix(samples, n_state, Mask_sample)

global e
e = 1e-15;

if nargin < 3
    Mask_sample = abs(samples) > 0;
end

if nargin < 2
    n_state = 2;     % 1: strong; 2: weak
end

alphas_state = zeros([1, n_state]);
mus_state = zeros([1, n_state]);
sigmas_state = zeros([1, n_state]);
N = sum(Mask_sample(:));


%% Initial guess
% thresh_arb = 0.1*max(abs(samples(:).*Mask_sample(:)))
% Mask_strong = ((samples.*Mask_sample) > thresh_arb);
% Mask_weak = Mask_sample & (~Mask_strong);
% alphas_state(1) = sum(Mask_strong(:))/N;
% alphas_state(2) = 1 - alphas_state(1);
% 
% mus_state(1) = mean(samples(Mask_strong));
% mus_state(2) = mean(samples(Mask_weak));
% 
% sigmas_state(1) = std(samples(Mask_strong))^2;
% sigmas_state(2) = std(samples(Mask_weak))^2;

% use {1~n-1}/n percentile to initialize
prctiles_state = zeros([1, n_state+1]);
prctiles_state(1) = -inf;
for i = 1:n_state-1
    prctiles_state(i+1) = prctile(samples(Mask_sample),i/n_state*100);
end
prctiles_state(end) = inf;
for i = 1:n_state
    Mask_tmp = Mask_sample & ...
               (samples.*Mask_sample) >= prctiles_state(i) & ...
               (samples.*Mask_sample) < prctiles_state(i+1);
    alphas_state(i) = sum(Mask_tmp(:))/N;
    mus_state(i) = mean(samples(Mask_tmp));
    sigmas_state(i) = std(samples(Mask_tmp))^2;
end


fprintf('Initial guess:\n');
show_val(n_state, alphas_state, mus_state, sigmas_state)


%% EM iteration

for i = 1:1000
   
    % Update alpha
    alphas_state_old = alphas_state;
    p_tmp = zeros([N, n_state]);
    for y = 1:n_state
        p_tmp(:, y) = calc_y_cond_x(y, n_state, samples(Mask_sample), alphas_state, mus_state, sigmas_state);
    end
    alphas_state = sum(p_tmp, 1)/N;
    
    % Update mu
    x_tmp = repmat(samples(Mask_sample), [1,n_state]);
    mus_state = sum(p_tmp.*x_tmp, 1)./sum(p_tmp, 1);
    
    % Update sigma
    var_tmp = (x_tmp - repmat(mus_state, [N,1])).^2;
    sigmas_state = sum(p_tmp.*var_tmp, 1)./sum(p_tmp, 1);
    
    % Output
    fprintf('Iter: %d\n', i);
    show_val(n_state, alphas_state, mus_state, sigmas_state)
    
    % Check if update is small enough
    if norm(alphas_state - alphas_state_old) < 1e-3
        break
    end
end



%% Find the optimal state for each voxel
p_tmp = zeros([N, n_state]);
for y = 1:n_state
    p_tmp(:, y) = calc_y_cond_x(y, n_state, samples(Mask_sample), alphas_state, mus_state, sigmas_state);
end
[dummy, idxs] = max(p_tmp, [], 2);

states = zeros([prod(size(samples)), 1]);
states(Mask_sample) = idxs;
states(~Mask_sample) = 1;     % Background is assigned to the first state
states = reshape(states, size(samples));

prob_states = zeros([size(samples), n_state]);
for y = 1:n_state
    tmp = zeros([prod(size(samples)), 1]);
    tmp(Mask_sample) = p_tmp(:, y);
    tmp(~Mask_sample) = (y==1);     % Background is assigned to the first state
    prob_states(:,:,:, y) = reshape(tmp, size(samples));
end


end


function p_x_cond_y = calc_x_cond_y(x, mu, sigma)

    global e

    % Gaussian
    p_x_cond_y = 1/(sqrt(2*pi)*sqrt(sigma))*exp(-1/2*(x-mu).^2/sigma);
    
    p_x_cond_y = max(p_x_cond_y, e);

end


function p_y_cond_x = calc_y_cond_x(y, n_state, x, alphas_state, mus_state, sigmas_state)

    % y is the state idx (1,2,...)
    
    p_x = 0;

    for i = 1:n_state
        p_x = p_x + alphas_state(i)*calc_x_cond_y(x, mus_state(i), sigmas_state(i));
    end
    
    p_y_cond_x = alphas_state(y)*calc_x_cond_y(x, mus_state(y), sigmas_state(y)) ./ p_x;

end


function show_val(n_state, alphas_state, mus_state, sigmas_state)
    
    fprintf('\t State\t alpha\t mu\t\t\t r(sigma)\n')
    for y = 1:n_state
        fprintf('\t\t %d\t %1.3f\t %1.1e\t %1.1e\n', y, alphas_state(y), mus_state(y), sqrt(sigmas_state(y)));
    end

end