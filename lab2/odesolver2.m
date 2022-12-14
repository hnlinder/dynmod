% fm = 1/2 *(-4 *a *alpha + alpha - 1)
% fa =  1/2 *(alpha *(-4* a - 4* m + 5) - 3)
% ga =1/2 *(alpha - 4 *alpha* m - 1)
% gm = 1/2 *(alpha* (-4* a - 4 *m + 5) - 3)
% 
% 
% 
% lamda = (fa - gm)/2 + sqrt((fa+gm)^2/4 - fa*gm + fm*ga);
% 
% 
% lamda_subbed = ((1/2 *(alpha *(-4* a - 4* m + 5) - 3)) - (1/2 *(alpha* (-4* a - 4 *m + 5) - 3)))/2 + sqrt(((1/2 *(alpha *(-4* a - 4* m + 5) - 3))+(1/2 *(alpha* (-4* a - 4 *m + 5) - 3)))^2/4 - (1/2 *(alpha *(-4* a - 4* m + 5) - 3))*(1/2 *(alpha* (-4* a - 4 *m + 5) - 3)) + (1/2 *(-4 *a *alpha + alpha - 1))*(1/2 *(alpha - 4 *alpha* m - 1)));




close all
clear
% lamda_subbed_a_to_m =@(alpha,m) ((1/2 *(alpha *(-4* m - 4* m + 5) - 3)) + (1/2 *(alpha* (-4* m - 4 *m + 5) - 3)))/2 + sqrt(((1/2 *(alpha *(-4* m - 4* m + 5) - 3))+(1/2 *(alpha* (-4* m - 4 *m + 5) - 3)))^2/4 - (1/2 *(alpha *(-4* m - 4* m + 5) - 3))*(1/2 *(alpha* (-4* m - 4 *m + 5) - 3)) + (1/2 *(-4 *m *alpha + alpha - 1))*(1/2 *(alpha - 4 *alpha* m - 1)));
m_of_alpha = @(alpha) (sqrt(3*alpha.^2 - 6*alpha + 4) + 3*alpha - 2)/(6*alpha);
lamda_subbed_a_to_m =@(alpha,m) ((1/2 *(alpha *(-4* m - 4* m + 5) - 3)) + (1/2 *(alpha* (-4* m - 4 *m + 5) - 3)))/2 + sqrt(((1/2 *(alpha *(-4* m - 4* m + 5) - 3))+(1/2 *(alpha* (-4* m - 4 *m + 5) - 3)))^2/4 - (1/2 *(alpha *(-4* m - 4* m + 5) - 3))*(1/2 *(alpha* (-4* m - 4 *m + 5) - 3)) + (1/2 *(-4 *m *alpha + alpha - 1))*(1/2 *(alpha - 4 *alpha* m - 1)));
alpha_of_F = @(F) F/(1+F);




F = linspace(0,10,1000000);
eigen = zeros([1,length(F)]);
for i=1:length(F)
    alpha = alpha_of_F(F(i));
    m = m_of_alpha(alpha);
    eigen(i) = lamda_subbed_a_to_m(alpha,m);
end
plot(F,eigen)
grid on 
posF = F(find(eigen >0));
posF(1)
%%
% fm = 1/2 *(-2 *a* alpha + alpha - 1) %yep
% fa = 1/2*(alpha *(-4 *a - 2* m + 5) - 3)%yep
% ga = 1/2 *(alpha - 2 *alpha*m - 1) %yep
% gm = 1/2 *(alpha*(-2*a - 4*m + 5) - 3)
% 
% 
% 
% lamda = (fa - gm)/2 + sqrt((fa+gm)^2/4 - fa*gm + fm*ga);
% 
% 
% lamda_subbed = ((1/2*(alpha *(-4 *a - 2* m + 5) - 3)) - (1/2 *(alpha*(-2*a - 4*m + 5) - 3)))/2 + sqrt(((1/2*(alpha *(-4 *a - 2* m + 5) - 3))+(1/2 *(alpha*(-2*a - 4*m + 5) - 3)))^2/4 - (1/2*(alpha *(-4 *a - 2* m + 5) - 3))*(1/2 *(alpha*(-2*a - 4*m + 5) - 3)) + (1/2 *(-2 *a* alpha + alpha - 1))*(1/2 *(alpha - 2 *alpha*m - 1)));
% lamda_subbed_a_to_m =@(alpha,m) ((1/2*(alpha *(-4 *m - 2* m + 5) - 3)) - (1/2 *(alpha*(-2*m - 4*m + 5) - 3)))/2 + sqrt(((1/2*(alpha *(-4 *m - 2* m + 5) - 3))+(1/2 *(alpha*(-2*m - 4*m + 5) - 3)))^2/4 - (1/2*(alpha *(-4 *m - 2* m + 5) - 3))*(1/2 *(alpha*(-2*m - 4*m + 5) - 3)) + (1/2 *(-2 *m* alpha + alpha - 1))*(1/2 *(alpha - 2 *alpha*m - 1)));

close all
clear

m_of_alpha =@(alpha) (sqrt(5*alpha.^2 - 8*alpha + 4) + 3 *alpha - 2)/(4 *alpha);
lamda_subbed_a_to_m =@(alpha,m) ((1/2*(alpha *(-4 *m - 2* m + 5) - 3)) - (1/2 *(alpha*(-2*m - 4*m + 5) - 3)))/2 + sqrt(((1/2*(alpha *(-4 *m - 2* m + 5) - 3))+(1/2 *(alpha*(-2*m - 4*m + 5) - 3)))^2/4 - (1/2*(alpha *(-4 *m - 2* m + 5) - 3))*(1/2 *(alpha*(-2*m - 4*m + 5) - 3)) + (1/2 *(-2 *m* alpha + alpha - 1))*(1/2 *(alpha - 2 *alpha*m - 1)));
alpha_of_F = @(F) F/(1+F);




F = linspace(0,100,1000000);
eigen = zeros([1,length(F)]);
for i=1:length(F)
    alpha = alpha_of_F(F(i));
    m = m_of_alpha(alpha);
    eigen(i) = lamda_subbed_a_to_m(alpha,m);
end
plot(F,eigen)
posF = F(find(eigen >0));
posF(1)




