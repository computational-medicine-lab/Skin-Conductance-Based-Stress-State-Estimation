function certainty = get_certainty_curve(vK, mu, xK, chance_prob)

    p = (1e-4:1e-4:1);
    [~, i] = min(abs(p - chance_prob));
    certainty = zeros(1, length(vK));
  
    for j = 1:length(vK)
        fp = cumtrapz(p, 1 ./ (sqrt(2 * pi * vK(j)) * p .* (1 - p)) .* ...
            exp(((-1) / (2 * vK(j)))* (log(p ./ ((1 - p) * exp(mu))) - xK(j)) .^ 2));
        certainty(1, j) = 1 - fp(i);
    end
end