function [lcl, ucl] = get_conf_lims(vkj, mu, xkj)

    p = (1e-4:1e-4:1);
  
    fp = cumtrapz(p, 1 ./ (sqrt(2 * pi * vkj) * p .* (1 - p)) .* ...
        exp(((-1) / (2 * vkj))* (log(p ./ ((1 - p) * exp(mu))) - xkj) .^ 2));
    
    n = find(fp <= 0.95);
    m = find(fp < 0.05);
    
    ucl = p(n(end));
    lcl = p(m(end));
end



