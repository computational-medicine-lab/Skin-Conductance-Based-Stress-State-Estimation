
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% FILTER-SMOOTHER FOR BINARY DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If you are using this code please cite the following paper

% Wickramasuriya DS, Amin MR and Faghih RT (2019) 
% Skin Conductance as a Viable Alternative for Closing the Deep Brain Stimulation Loop in Neuropsychiatric Disorders. 
% Front. Neurosci. 13:780. doi: 10.3389/fnins.2019.00780

close all;
clear;
clc;
fs = 2;
fsu = 4;
fss = 2;
bw = fsu;
subjects = [1, 5, 8, 9, 12, 16];

all_chance_probs = zeros(1, length(subjects));
k = 1;

% load('s_iir_filt_cvx_eda.mat');
mini_emo_start_i = 4;
cog_stress_start_i = 5;
relax_start_i = 6;
emo_stress_start_i = 7;
final_relax_start_i = 8;

min_peak_distance = 1 * fsu;
min_peak_height = 0.35;
min_peak_prominence = 0.2;

end_ = 5*60*fsu;
addpath('Brewer_Map_Color_library');

% my colors
cg = [163 255 101]/255;
cb = [232 190 102]/255;
cr = [255 73 72]/255;
cv = [126 97 232]/255;
cc = [68 255 239]/255;

for i = subjects
    data = load(['deconvolution_results\result_s',num2str(i),'.mat']);
    
    u = (data.uj_whole(1:end-end_)'>0);
    u = u(:)';
    
    n = [];
    j = 1;

    while (j * bw) <= length(u)
        if sum(u((1 + (j - 1) * bw):(j * bw))) > 0
            n = [n 1];
        else
            n = [n 0];
        end
        
        j = j + 1;
    end
    
    all_chance_probs(k) = sum(n) / length(n);
    k = k + 1;
end

for i = subjects
    sub = i;
    data = load(['deconvolution_results\result_s',num2str(i),'.mat']);

    u = (data.uj_whole(1:end-end_)'>0); u = u(:)';

    n = [];
    j = 1;

    while (j * bw) <= length(u)
        if sum(u((1 + (j - 1) * bw):(j * bw))) > 0
            n = [n 1];
        else
            n = [n 0];
        end
        j = j + 1;
    end

    permit_bias = false;
    x0_guess = 0;
    veps_guess = 0.005;   % eps_k variance estimate
    chance_prob = mean(all_chance_probs);
    chance_prob = sum(n)/length(n);
    % chance_prob = all_chance_probs(subjects == i);
    tol = 1e-8;
    M = 20000; % num itr

    x0 = zeros(1, M + 1);
    veps = zeros(1, M + 1);
    N = length(n);
    mu = log(chance_prob / (1 - chance_prob));

    veps(1) = veps_guess;
    x0(1) = x0_guess;

    for m = 1:M

        xk = zeros(N, N); vk = zeros(N, N);    
        xK = zeros(1, N); vK = zeros(1, N); 

        cvK = zeros(N, N); WK = zeros(N, N);    
        AK = zeros(1, N); 

        for k = 1:N
            if (k == 1)                 
                if (m == 1)
                    xk(1, 1) = get_posterior_mode(mu, veps_guess + veps(m), n(k), x0(m));
                    vk(1, 1) = 1 / ((1 / (veps_guess + veps(m))) + (exp(mu + xk(1, 1)) / (1 + exp(mu + xk(1, 1)))) ...
                        * (1 / (1 + exp(mu + xk(1, 1)))));
                else
                    xk(1, 1) = get_posterior_mode(mu, veps(m - 1) + veps(m), n(k), x0(m));
                    vk(1, 1) = 1 / ((1 / (veps(m - 1) + veps(m))) + (exp(mu + xk(1, 1)) / (1 + exp(mu + xk(1, 1)))) ...
                        * (1 / (1 + exp(mu + xk(1, 1)))));  
                end
            else 
                xk(k, k - 1) = xk(k - 1, k - 1);  
                vk(k, k - 1) = vk(k - 1, k - 1) + veps(m);

                xk(k, k) = get_posterior_mode(mu, vk(k, k - 1), n(k), xk(k, k - 1));
                vk(k, k) = 1 / ((1 / vk(k, k - 1)) + (exp(mu + xk(k, k)) / (1 + exp(mu + xk(k, k)))) ...
                    * (1 / (1 + exp(mu + xk(k, k)))));        
            end
        end

        xK(1, N) = xk(N, N);
        vK(1, N) = vk(N, N);

        for k = (N - 1):(-1):1
            AK(k) = vk(k, k) / vk(k + 1, k);
            xK(k) = xk(k, k) + AK(k) * (xK(k + 1) - xk(k + 1, k));
            vK(k) = vk(k, k) + AK(k) * AK(k) * (vK(k + 1) - vk(k + 1, k));
        end

        for k = 1:N
            cvK(k, k) = vK(k);

            if (k > 1)
                cvK(k - 1, k) = AK(k - 1) * cvK(k, k);
            end

            WK(k, k) = vK(k) + (xK(k) ^ 2); 

            if (k > 1)  
                WK(k - 1, k) = cvK(k - 1, k) + xK(k) * xK(k - 1);
            end        
        end

        if (m < M) 
            if (~permit_bias)            
                x0(m + 1) = 0.5 * xK(1);
                veps(m + 1) = 1 / (N + 1) * (2 * (sum(diag(WK)) - WK(1, 1)) ...
                    - 2 * sum(diag(WK, 1)) + (1.5 * xK(1) ^ 2 + 2 * vK(1)) - WK(N, N));
            else
                x0(m + 1) = xK(1);
                veps(m + 1) = (1 / (N + 1)) * (2 * (sum(diag(WK)) - WK(1, 1)) ...
                    - 2 * sum(diag(WK, 1)) + (1.5 * xK(1) ^ 2 + 2 * vK(1)) - WK(N, N));
            end

    %         fprintf('\nx(%d) = %.18f\nv(%d) = %.18f\n', (m + 1), x0(1, m + 1), ...
    %             (m + 1), veps(1, m + 1));

            if (abs(veps(m + 1) - veps(m)) < tol) && (abs(x0(m + 1) - x0(m)) < tol) 
                fprintf('\n\nConvergence after %d steps x0 = %.10f v_eps = %.10f\n\n', m, x0(m), veps(m));
                break;
            end
        end

    end

    pk = zeros(1, N); pK = zeros(1, N);
    lclk = zeros(1, N); uclk = zeros(1, N);
    lclK = zeros(1, N); uclK = zeros(1, N);
    fp = zeros(N, 1e4 + 1);

    fprintf('\n\nPlotting\n');

    for k = 1:N
        pk(k) = get_fp_mode(vk(k, k), mu, xk(k, k));
        [pK(k), fp(k, :), y] = get_fp_mode(vK(k), mu, xK(k));

        [lclK(k), uclK(k)] = get_conf_lims(vK(k), mu, xK(k));
        [lclk(k), uclk(k)] = get_conf_lims(vk(k, k), mu, xk(k, k));   
        fprintf('.');

        if mod(k, 100) == 0
            fprintf('\n');
        end
    end
    fprintf('\n\n');
    save(['stress_estimation_results\result_stress_',num2str(sub),'.mat']);
end
