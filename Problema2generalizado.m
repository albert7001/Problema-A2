% código para el modelo 1D teórico con N>>, mismo formato que el exacto
format long

syms eps mu B b k x

N = 4;

Z = exp(N*b*eps)*(cosh(b*mu*B)+sqrt(cosh(b*mu*B)^2-2*exp(-2*b*eps)*sinh(2*b*eps)))^N;

E_med = -diff(log(Z),b);

F = -1/b*log(Z);

S = k*(log(Z)+b*E_med);

M = 1/b*diff(log(Z),B);

%

eps_s = 20;
mu_s = 1;
B_s = 1;
b_s = 0.02;
k_s = 0.1;

Z_b = subs(Z,{eps,mu,B},{eps_s,mu_s,B_s});
Z_s = double(subs(Z_b,{b},{b_s}));

E_med_b = subs(E_med,{eps,mu,B},{eps_s,mu_s,B_s})+N*eps_s+N*mu_s*B_s;
E_med_s = double(subs(E_med_b,{b},{b_s}));

M_B1 = subs(M,{eps,mu,b},{eps_s,mu_s,b_s});
M_B2 = subs(M,{eps,mu,b},{eps_s,mu_s,0.04});
M_s = double(subs((M_B1),{B},{B_s}));

S_b = subs(S,{eps,mu,B,k},{eps_s,mu_s,B_s,k_s});
S_s = double(subs(S_b,{b},{b_s}));
Smax = double(subs(S_b,{b},{0}));

m = (1-(sinh(2/x))^-4)^(1/8);


figure
fplot(Z_b,[0,0.4],'LineWidth',1.5)
xlabel('$\beta$','Interpreter','latex')
ylabel('$Z$','Interpreter','latex')
set(gcf,'Color','w')

figure % las unidades de magnetización en 1D son A*m^2, momento dipolar magnético por unidad de espín
fplot(M_B1,[-100,100],'LineWidth',1.5)
hold on
fplot(M_B2,[-100,100],'LineWidth',1.5)
xticks(0)
yticks([-N 0 N])
yticklabels({'$-N\mu_m$','0','$N\mu_m$'})
xlabel('$B$','Interpreter','latex')
ylabel('$M$','Interpreter','latex')
set(gca, 'TickLabelInterpreter', 'latex')
set(gcf,'Color','w')
legend('T_1','T_2<T_1')

figure
fplot(S_b/Smax,[0,1],'LineWidth',1.5)
xlabel('$\beta$','Interpreter','latex')
ylabel('$S/S_{max}$','Interpreter','latex')
xticks(0:0.25:1)
yticks([0 1])
set(gcf,'Color','w')

figure
fplot(m,[0,2.269],'LineWidth',1.5)
hold on
fplot(-m,[0,2.269],'LineWidth',1.5,'Color',[0 0.4470 0.7410])
line([2.269 2.269], [-0.377 0.377],'LineWidth',1.5,'Color',[0 0.4470 0.7410])
line([2.269 4], [0 0],'LineWidth',1.5,'Color',[0 0.4470 0.7410])
xlabel('$T$','Interpreter','latex')
ylabel('$m$','Interpreter','latex')
xticks([0 2.269])
xticklabels({'0','$T_c$'})
yticks([-1 0 1])
yticklabels({'$-\mu_m$','0','$\mu_m$'})
xlim([0 4])
set(gca, 'TickLabelInterpreter', 'latex')
set(gcf,'Color','w')
saveas(gcf,'mag.png')