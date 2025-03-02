format long

syms eps mu B b k

% número de espines (este código permite introducir un número arbitrario de
% espines)
N = 4;

% las filas son los estados y las columnas el valor del espín correspondiente
A = zeros(2^N,N);

% bucle que extrae los números en binario en una matriz
for i = 0: 2^N-1
    A(i+1,:) = bitget(i,N:-1:1);
end

A(A==0) = -1; % cambio de 0 a -1

% definición de la energía (se suma cada contribución en un bucle)
E_ne = 0;
for i = 1: N
    if i == N
        E_ne = E_ne-eps*(A(:,i).*A(:,1))-mu*B*A(:,i);
    else
    E_ne = E_ne-eps*(A(:,i).*A(:,i+1))-mu*B*A(:,i);
    end
end

% reescalado de la energía (fundamental E=0)
E = E_ne+N*eps+N*mu*B;

% función de partición
Z = sum(exp(-b*E));

% energía media
E_med = -diff(log(Z),b);

% función de Helmholtz
F = -1/b*log(Z);

% entropía
S = k*(log(Z)+b*E_med);

% magnetización
M = 1/b*diff(log(sum(exp(-b*E_ne))),B);

M_teo = (N*mu*sinh(b*mu*B))/sqrt(cosh(b*mu*B)^2-2*exp(-2*b*eps)*sinh(2*b*eps));


% constantes ajustables
eps_s = 20;
mu_s = 1;
B_s = 1;
b_s = 0.02;
k_s = 0.1;

% sustitución por las constantes y funciones de una variable

E_s = double(subs(E,{eps,mu,B},{eps_s,mu_s,B_s}));

Z_b = subs(Z,{eps,mu,B},{eps_s,mu_s,B_s});
Z_s = double(subs(Z_b,{b},{b_s}));

E_med_s = double(subs(E_med,{eps,mu,B,b},{eps_s,mu_s,B_s,b_s}));

F_b = subs(F,{eps,mu,B},{eps_s,mu_s,B_s});
F_s = double(subs(F_b,{b},{b_s}));

S_b = subs(S,{eps,mu,B,k},{eps_s,mu_s,B_s,k_s});
S_s = double(subs(S_b,{b},{b_s}));
Smax = double(subs(S_b,{b},{0}));

M_B1 = subs(M,{eps,mu,b},{eps_s,mu_s,b_s});
M_B2 = subs(M,{eps,mu,b},{eps_s,mu_s,0.04});
M_s = double(subs(M_B1,{B},{B_s}));
M_b_teo = subs(M_teo,{eps,mu,b,k},{eps_s,mu_s,b_s,k_s});
M_teo_s = double(subs(M_b_teo,{B},{B_s}));

% energías
figure
stem(E_s,'filled','LineWidth',1.5);
xlabel('Microestado')
ylabel('$E$/ au','Interpreter','latex')
set(gcf,'Color','w')

% función de partición
figure
fplot(Z_b,[0,0.4],'LineWidth',1.5)
xlabel('$\beta$','Interpreter','latex')
ylabel('$Z$','Interpreter','latex')
set(gcf,'Color','w')

% entropía
figure
fplot(S_b/Smax,[0,1],'LineWidth',1.5)
xlabel('$\beta$','Interpreter','latex')
ylabel('$S/S_{max}$','Interpreter','latex')
xticks(0:0.25:1)
yticks([0 1])
set(gcf,'Color','w')

% magnetización
figure
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