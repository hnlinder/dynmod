
%% solving odes in lab2
clear
close all
F = 1;
alpha = F/(1+F);
L = 60;
tmax = 10000;
t_span = [0,L*tmax];

Y0 = [29, 29];

opts = odeset('RelTol',1e-6,'AbsTol',1e-4);

[t,y] = ode45(@(t,Y) odesolver_func(t,Y,alpha) , t_span , Y0, opts);
% [t113,y113] = ode113(@(t,Y) oregonator_func(t,Y,alpha,beta,gamma) , t_span , Y0, opts);
% [t23s,y23s] = ode23s(@(t,Y) oregonator_func(t,Y,alpha,beta,gamma) , t_span , Y0, opts);
% [t15s,y15s] = ode15s(@(t,Y) oregonator_func(t,Y,alpha,beta,gamma) , t_span , Y0, opts);
%%
close all 
figure 
plot(t,y)
   % phase plot
   figure%(2)
   hold on
   plot(y(:,1), y(:,2))
   xlabel("a")
   ylabel("m")
%    title(sprintf("rsat = %d, rgrass = %d", r_sat, r_grass))
%    legend
% %%
% tstep = [t(1:end-1) - t(2:end)];
% close all
% figure%(1)
% hold on
% 
% % yyaxis left
% plot(t,y45)
% %   plot(t15s,y15s)
% %   plot(t113,y113)
% %   plot(t23s,y23s)
% % yyaxis right
% % plot(t(1:end-1), tstep)
% xlabel("t")
% ylabel("y")
% %  title(sprintf("rsat = %d, rgrass = %d", r_sat, r_grass))
% legend("45", "15s", "113", "23s")
% %%
% spacing = 1;
% xmax = 60;
% [x,y] = meshgrid(-0:spacing:xmax,-0:spacing:xmax);
% % dydt = brusselator_func(t,[x;y], a, b);
% % scalefactor = .01;
% % dx = dydt(1:length(dydt(1,:)),:);
% % dy = dydt(length(dydt(1,:))+1:end,:);
% % ndx = dx./(sqrt(dx.^2 + dy.^2)) * scalefactor;
% % ndy = dy./(sqrt(dx.^2 + dy.^2)) * scalefactor;
% 
% % phase plot
% figure%(2)
% hold on
% plot3(y45(:,1), y45(:,2), y45(:,3))
% % quiver(x,y,ndx,ndy,0)
% xlabel("x")
% ylabel("y")
% %  title(sprintf("rsat = %d, rgrass = %d", r_sat, r_grass))
% legend
% %%
% % figure
% hold on 
% x = x(1,:);
% plot(x,b/a./x,"LineWidth",2)
% plot(x,((b + 1).*x - 1)./(a*x.^2), "LineWidth",2)
% xlim([-0,3])
% ylim([-0,3])
% legend("Spiral","$(\dot x, \dot y )$","$x = \frac{b}{ay}$","$y = \frac{((b + 1)x - 1)}{(ax^2)}$","Interpreter","Latex")
%   %[t,y15s] = ode15s(@(t,Y) brusselator_func(t,Y,a,b) , t_span , Y0, opts);

