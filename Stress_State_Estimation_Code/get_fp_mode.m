function [p_mode, fp, p] = get_fp_mode(vkj, mu, xkj)

    p = (0:0.0001:1);
    syms p_sym;
    
    fp = 1 ./ (sqrt(2 * pi * vkj) * p .* (1 - p)) .* ...
        exp(((-1) / (2 * vkj))* (log(p ./ ((1 - p) * exp(mu))) - xkj) .^ 2);

%     p_mode = double(vpasolve(((2 * p_sym - 1) - 1 / vkj * ...
%         (-xkj + log(p_sym) - log(1 - p_sym) - mu)) == 0, p_sym, p(fp == max(fp))));
    
    [~, i] = max(fp);
    p_mode = p(i);

    %fprintf('\np = %.10f\np = %.10f\n', p_mode, p_mode2);
end

