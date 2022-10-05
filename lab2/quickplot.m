m = @(alpha) (sqrt(3*alpha.^2 - 6*alpha + 4) - 3*alpha + 2)./(6*alpha)


close all
m(10)
fplot(m, [0,1])

x = linspace(0,3,10000);

x(find(m(x)<0))(1)
