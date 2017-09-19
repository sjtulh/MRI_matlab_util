% Non-linear Conjugate CG with Back-Tracking Line Search
%   E = lambda*||Gx||_1 + (|| [DCF].^0.5 .* (Fx - b)||_2^2)

function [x, res_gE_rel, k, x_hist, res_gE_rel_hist, cost_reg_hist, cost_fidelity_hist, norm_res_gE_fidelity_hist] = cgsolve_mod(f_reg, f_fidelity, b_tilde, f_E_reg, f_E_fidelity, tol_grad, maxiter, verbose, x0, tol_gE_fidelity, flag_skip_LS)

tic

if nargin < 11
    flag_skip_LS = 0;
end

matrix_size = size(x0);
x = x0;

k = 0;

if ~flag_skip_LS
    E_reg = f_E_reg(x);
    E_fidelity = f_E_fidelity(x);
else
    E_reg = 0;
    E_fidelity = 0;
end
E_k =  E_reg + E_fidelity;
gE_fidelity = f_fidelity(x);
gE_k = f_reg(x,x) + gE_fidelity - b_tilde;
p_k = -gE_k;
gE_k = gE_k(:);

gE_0 = gE_k;
res_gE_rel = (gE_k'*gE_k)/(gE_0'*gE_0);

% Log
x_hist(:,:,:,k+1) = x;
cost_reg_hist(k+1) = E_reg;
cost_fidelity_hist(k+1) = E_fidelity;
res_gE_rel_hist(k+1) = res_gE_rel;
res_gE_fidelity = gE_fidelity - b_tilde;
res_gE_fidelity_rel = res_gE_fidelity(:)'*res_gE_fidelity(:)/(b_tilde(:)'*b_tilde(:));
norm_res_gE_fidelity_hist(k+1) = res_gE_fidelity(:)'*res_gE_fidelity(:);

toc

while res_gE_rel > tol_grad && k < maxiter && res_gE_fidelity_rel > tol_gE_fidelity
    
    % Parameters for Back-Tracking Line Search
    s = 0.05;
    beta = 0.6;

    % Initialize alpha_k
%     alpha_k = 1;
    % Use alpha from CG
    q_fidelity = f_fidelity(p_k);
    q = f_reg(p_k,x) + q_fidelity; 
    q = q(:);
    alpha_k = (gE_k'*gE_k)/(p_k(:)'*q);
    alpha_k = real(alpha_k);
    
    
    if ~flag_skip_LS
        % Back-Tracking Line Search
        E_reg = f_E_reg(x + alpha_k*p_k);
        E_fidelity = f_E_fidelity(x + alpha_k*p_k);
        E_kp1 = E_reg + E_fidelity;
        while E_kp1 > (E_k + s*alpha_k*real(gE_k'*p_k(:)))
            disp(sprintf('alpha_k = %8.3e failed, shrink it by %f', alpha_k, beta));
            alpha_k = alpha_k*beta;
            E_reg = f_E_reg(x + alpha_k*p_k);
            E_fidelity = f_E_fidelity(x + alpha_k*p_k);
            E_kp1 = E_reg + E_fidelity;   
        end
    else
        % Skip Line Search (for speed)
        E_reg = 0;
        E_fidelity = 0;
    end
    
    % Update x
    %   x_k+1 = x_k + alpha_k*p_k;
    x = x + alpha_k*p_k;
    
    
    % Evaluate gE_k+1
    if (mod(k+1,50) == 0)
        gE_fidelity = f_fidelity(x);
        gE_kp1 = f_reg(x,x) + gE_fidelity - b_tilde;
        gE_kp1 = gE_kp1(:);
        res_gE_rel = (gE_kp1'*gE_kp1)/(gE_0'*gE_0);
    else
        gE_fidelity = gE_fidelity + alpha_k*q_fidelity;
        gE_kp1 = f_reg(x,x) + gE_fidelity - b_tilde;
        gE_kp1 = gE_kp1(:);
        res_gE_rel = (gE_kp1'*gE_kp1)/(gE_0'*gE_0);
    end
    
    % Update p_k
    beta_FR = (gE_kp1'*gE_kp1)/(gE_k'*gE_k);
    p_kp1 = -reshape(gE_kp1, matrix_size) + beta_FR * p_k;
    
    % Update k
    k = k+1;
    p_k = p_kp1;
    gE_k = gE_kp1;
      
    % Log
    x_hist(:,:,:,k+1) = x;
    cost_reg_hist(k+1) = E_reg;
    cost_fidelity_hist(k+1) = E_fidelity;
    res_gE_rel_hist(k+1) = res_gE_rel;
    res_gE_fidelity = gE_fidelity - b_tilde;
    res_gE_fidelity_rel = res_gE_fidelity(:)'*res_gE_fidelity(:)/(b_tilde(:)'*b_tilde(:));
    norm_res_gE_fidelity_hist(k+1) = res_gE_fidelity(:)'*res_gE_fidelity(:);   
    
    % Display
    if verbose
        disp(sprintf('cg: Iter = %d, Residual = %8.3e, gE fidelity = %8.3e', k, sqrt(res_gE_rel), sqrt(res_gE_fidelity_rel)));
        disp(sprintf('              cost reg = %8.3e, cost fidelity = %8.3e', E_reg, E_fidelity));
    end
    
    toc
    
end
    









