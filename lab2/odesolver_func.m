function dydt = odesolver_func(t, x_vec,alpha)
  a = x_vec(1, :);
  m = x_vec(2, :);
  dydt = [
          alpha*(a.*(1-a-m) - a.*m) + (1 - alpha).*((1-a-m)/2 - a);
          alpha*(m.*(1-a-m) - a.*m) + (1 - alpha).*((1-a-m)/2 - m)
          ];
end