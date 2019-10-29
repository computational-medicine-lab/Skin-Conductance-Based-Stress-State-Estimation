function [sol] = get_posterior_mode(mu, vkj, nk, xkj)
  
    it(1) = xkj + vkj * (nk - exp(mu + xkj) / (1 + exp(mu + xkj)));

    for i = 1:40    
       g(i)     = xkj + vkj * (nk - exp(mu + it(i)) / (1 + exp(mu + it(i)))) - it(i);
       gprime(i)= (-1) * vkj * exp(mu + it(i)) / (1 + exp(mu + it(i)))^2 - 1;
       it(i+1)  = it(i) - g(i)/gprime(i);

       sol = it(i+1);
       if abs(sol - it(i))<1e-14    
          return
       end
    end    

end

