%% Nova Flexao Normal Composta


% Limpando variaveis
clear all
close all
clc

% Especificacoes concreto
fck = 60; %MPa
gamma_c = 1.4;
sigma_cd = 0.85*fck/gamma_c; %MPa
gamma_conc = 2500; %kg/m3

% Especificacoes aco
f_yk = 500; %Mpa
gamma_s = 1.15;
E_s = 210; %GPa
f_yd = f_yk/gamma_s; %MPa
epsilon_yd = f_yd/E_s; % por mil
gamma_aco = 7850; %kg/m3


% Especificacoes de design
c = 0.04; % m
b = 0.25; % m
h = 0.70; % m
d = 0.05; % m
nc = 3;
nb = [2 2 2];
phi = 0.020; % m
y_s = [-0.30, -0.22, 0.30];

% Esforcos solicitantes

N_d = 4.930*1.4; % MN # Majorar esforcos
M_d = 0*1.4; % MN*m # Majorar esforcos

% Inicializando variaveis

epsilon_c2 = 0; % por mil
epsilon_cu = 0; % por mil
x_lim = 0;
n = 0;

% Tolerancias

tol_J = 1e-9;
tol_k = 1e-9;
tol_f = 1e-9;
i = 0;
it_max = 1e6;

% Especificacoes do problema
if fck <= 50
    x_lim = 0.45*d;
    n = 2;
    epsilon_c2 = 2;
    epsilon_cu = 3.5;
else
    if (50 < fck) && (fck <= 90)
        x_lim = 0.35*d;
        n = 1.4 + 23.4*power((90-fck)/100,4);
        epsilon_c2 = 2 + 0.085*power(fck-50,0.53);
        epsilon_cu = 2.6 + 35*power((90-fck)/100,4);
    end
end

%% Plotando os esforços resistentes

y_t=h/2;
y_b=-h/2;
% plot_epsilon_0_k(epsilon_c2,epsilon_cu,h,y_t,y_b,y_s)
% plot_N_r_M_r(epsilon_c2,epsilon_cu,sigma_cd,n,b,f_yd,epsilon_yd,h,y_t,y_b,y_s,phi,nb,tol_k)

%% Testando um valor especifico
epsilon_0 = 0.2734;
k = -7.7542;


N_s = Ns(y_s,nb,phi,f_yd,epsilon_yd,epsilon_0,k);
M_s = Ms(y_s,nb,phi,f_yd,epsilon_yd,epsilon_0,k);

N_c = Nc(epsilon_0, k, epsilon_c2, sigma_cd,n,y_b,y_t,b,h,tol_k);
M_c = Mc(epsilon_0, k, epsilon_c2, sigma_cd,n,y_b,y_t,b,h,tol_k);

N_r = double(N_c + N_s);
M_r = double(M_c + M_s);

epsilon_t = epsilon_0 + y_t*k;
epsilon_b =  epsilon_0 + y_b*k;


%% Metodo de Newton Raphson

epsilon_0 = 0;
k = 0;
epsilon_0_it = epsilon_0;
k_it = k;

N_s = Ns(y_s,nb,phi,f_yd,epsilon_yd,epsilon_0,k);
M_s = Ms(y_s,nb,phi,f_yd,epsilon_yd,epsilon_0,k);

N_c = Nc(epsilon_0, k, epsilon_c2, sigma_cd,n,y_b,y_t,b,h,tol_k);
M_c = Mc(epsilon_0, k, epsilon_c2, sigma_cd,n,y_b,y_t,b,h,tol_k);

N_r = double(N_c + N_s);
M_r = double(M_c + M_s);

% Norma euclidiana admensionalizada
f_ad = sqrt(power((N_d-N_r)/(sigma_cd*b*h),2)+power((M_d-M_r)/(sigma_cd*b*h*h),2));
while (abs(f_ad) > tol_f) && (i < it_max)
    i = i + 1;
    epsilon_t = epsilon_0 + k*y_t;
    epsilon_b = epsilon_0 + k*y_b;
    
    N_s = Ns(y_s,nb,phi,f_yd,epsilon_yd,epsilon_0,k);
    M_s = Ms(y_s,nb,phi,f_yd,epsilon_yd,epsilon_0,k);

    N_c = Nc(epsilon_0, k, epsilon_c2, sigma_cd,n,y_b,y_t,b,h,tol_k);
    M_c = Mc(epsilon_0, k, epsilon_c2, sigma_cd,n,y_b,y_t,b,h,tol_k);
    
    N_r = double(N_c + N_s);
    M_r = double(M_c + M_s);
 
    J_c = Jc(epsilon_t,epsilon_b, epsilon_0, epsilon_c2, sigma_cd,n,b,k,h,tol_k);
    J_s = Js(E_s,epsilon_yd,phi,y_s,nb,epsilon_0,k);
    J = J_c + J_s; % J = -[EA ES; ES EI]
    J_ad = J(1,1)*J(2,2)/(sigma_cd*b*h*sigma_cd*b*h*power(h,2)) - (J(1,2)*J(2,1))/power(sigma_cd*b*h*h,2);
    if abs(J_ad) > tol_J
        J_inv = power(J(1,1)*J(2,2)-J(1,2)*J(2,1),-1)*[J(2,2), -J(2,1); -J(1,2), J(1,1)];
        f = [N_d - N_r; M_d - M_r];
        new_it = [epsilon_0; k] - J_inv*f;
        epsilon_0 = new_it(1);
        epsilon_0_it = [epsilon_0_it, epsilon_0];
        k = new_it(2);
        k_it = [k_it, k];
        f_ad = sqrt(power((N_d-N_r)/(sigma_cd*b*h),2)+power((M_d-M_r)/(sigma_cd*b*h*h),2));
    else
        fprintf("Nao existe solucao!\n")
        break
    end
end
ELU(epsilon_0,k,y_b,y_t,y_s,epsilon_c2,epsilon_cu)

plot_epsilon_0_k(epsilon_c2,epsilon_cu,h,y_t,y_b,y_s,epsilon_0_it,k_it)
plot_N_r_M_r(epsilon_c2,epsilon_cu,sigma_cd,n,b,f_yd,epsilon_yd,h,y_t,y_b,y_s,phi,nb,tol_k)

%% Funcoes de plot

function plot_epsilon_0_k(epsilon_c2,epsilon_cu,h,y_t,y_b,y_s,epsilon_0_it,k_it)
c = epsilon_c2;
u = epsilon_cu;
A = [0;c];
B = [u/h;u*y_t/h];
C = [(u+10)/(y_t-min(y_s));u-y_t*((u+10)/(y_t-min(y_s)))];
D = [0;-10];
E = [(u+10)/(y_b-max(y_s));u-y_b*((u+10)/(y_b-max(y_s)))];
F = [-u/h;-u*y_b/h];
points = [A B C D E F A];
figure(1);
set(gca, 'FontSize', 12);
set(gca, 'TickLabelInterpreter', 'latex');
plot(points(1,:),points(2,:),'-b', 'LineWidth', 2);
% Plotar eixos
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';xlabel('$\kappa \; (1/m)$', Interpreter='latex')
ylabel('$\epsilon_0$', Interpreter='latex')
xlim(1.1*[min(points(1,:)) max(points(1,:))])
ylim(1.1*[min(points(2,:)) max(points(2,:))])
title('Regiao viavel para o dimensionamento', Interpreter='latex')

hold on
plot(k_it,epsilon_0_it,'-g')
plot(k_it(length(k_it)),epsilon_0_it(length(epsilon_0_it)),'-gx')
hold off
end

function plot_N_r_M_r(epsilon_c2,epsilon_cu,sigma_cd,n,b,f_yd,epsilon_yd,h,y_t,y_b,y_s,phi,nb,tol_k)
c = epsilon_c2;
u = epsilon_cu;
A = [0;c];
B = [u/h;u*y_t/h];
C = [(u+10)/(y_t-min(y_s));u-y_t*((u+10)/(y_t-min(y_s)))];
D = [0;-10];
E = [(u+10)/(y_b-max(y_s));u-y_b*((u+10)/(y_b-max(y_s)))];
F = [-u/h;-u*y_b/h];
points = [A B C D E F A];

num_points = 1000;
kappa_values = [];
epsilon_0_values = [];
for i=1:6
    kappa_values = [kappa_values, linspace(points(1,i), points(1,i+1), num_points)];
    epsilon_0_values = [epsilon_0_values,  linspace(points(2,i), points(2,i+1), num_points)];

end
N_r_values = zeros(1,length(kappa_values));
M_r_values = zeros(1,length(epsilon_0_values));

for j=1:length(kappa_values)
    N_c = Nc(epsilon_0_values(j), kappa_values(j), epsilon_c2, sigma_cd,n,y_b,y_t,b,h,tol_k);
    M_c = Mc(epsilon_0_values(j), kappa_values(j), epsilon_c2, sigma_cd,n,y_b,y_t,b,h,tol_k);
    N_s = Ns(y_s,nb,phi,f_yd,epsilon_yd,epsilon_0_values(j),kappa_values(j));
    M_s = Ms(y_s,nb,phi,f_yd,epsilon_yd,epsilon_0_values(j),kappa_values(j));
    N_r_values(j) = (N_c + N_s)*1e3;
    M_r_values(j) = (M_c + M_s)*1e3;
end

figure(2);
plot(N_r_values, M_r_values, '-r', 'LineWidth', 2);
% Plotar eixos
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';xlabel('$N_r \; (kN)$', 'Interpreter', 'latex');
ylabel('$M_r \; (kNm)$', 'Interpreter', 'latex');
xlim([min(N_r_values)*1.1, max(N_r_values)*1.1]);
ylim([min(M_r_values)*1.1, max(M_r_values)*1.1]);
title('Esforcos Resistentes', 'Interpreter', 'latex');
set(gca, 'FontSize', 12);
set(gca, 'TickLabelInterpreter', 'latex');

end

%% Funcoes potencial

% Nao usada
% function I_m = Im(epsilon, epsilon_c2, sigma_cd,n,m)
% syms csi
% if epsilon <= 0
%     I_m = 0;
% else
%     if (0 < epsilon) && (epsilon <= epsilon_c2)
%         I_m = int(sigma_cd*(1-power(1-csi/epsilon_c2,n))*power(csi,m),csi,0,epsilon);
%     else
%         if epsilon > epsilon_c2
%             I_m = int(sigma_cd*(1-power(1-csi/epsilon_c2,n))*power(csi,m),csi,0,epsilon_c2) + sigma_cd*int(power(csi,m),csi,epsilon_c2,epsilon);
%         end
%     end
% end
% end

function I_0 = I0(epsilon, epsilon_c2, sigma_cd,n)
    if epsilon < 0
        I_0 = 0;
    else
    if (0 <= epsilon) && (epsilon <= epsilon_c2)
        I_0 = sigma_cd*(epsilon - epsilon_c2*(1 - power(1 - epsilon/epsilon_c2,n + 1))/(n + 1));
    else
        I_0 = sigma_cd*(epsilon - epsilon_c2/(n + 1));
    end
    end
end

function I_1 = I1(epsilon, epsilon_c2, sigma_cd,n)
    if epsilon < 0
        I_1 = 0;
    else
    if (0 <= epsilon) && (epsilon <= epsilon_c2)
        I_1 = sigma_cd*(power(epsilon,2)/2 + power(epsilon_c2,2)*(-1*(1 - power(1 - epsilon/epsilon_c2,n + 1))/(n + 1) + 1*(1 - power(1 - epsilon/epsilon_c2,n + 2))/(n + 2)));

    else
        I_1 = sigma_cd*(power(epsilon,2)/2 - power(epsilon_c2,2)/((n + 1)*(n + 2)));
    end
    end
end

% Nao usada
% function DI_m = DIm(epsilon_t,epsilon_b, epsilon_c2, sigma_cd,n,m)
% DI_m = Im(epsilon_t, epsilon_c2, sigma_cd,n,m)-Im(epsilon_b, epsilon_c2, sigma_cd,n,m);
% end

function DI_0 = DI0(epsilon_t,epsilon_b, epsilon_c2, sigma_cd,n)
DI_0 = I0(epsilon_t, epsilon_c2, sigma_cd,n)-I0(epsilon_b, epsilon_c2, sigma_cd,n);
end


function DI_1 = DI1(epsilon_t,epsilon_b, epsilon_c2, sigma_cd,n)
DI_1 = I1(epsilon_t, epsilon_c2, sigma_cd,n)-I1(epsilon_b, epsilon_c2, sigma_cd,n);
end

function DJ_m = DJm(epsilon_t,epsilon_b, epsilon_0, epsilon_c2, sigma_cd,n,m)
DJ_m = Jm(epsilon_t, epsilon_0, epsilon_c2, sigma_cd,n,m)-Jm(epsilon_b, epsilon_0, epsilon_c2, sigma_cd,n,m);
end

function J_m = Jm(epsilon, epsilon_0, epsilon_c2,sigma_cd,n,m)
J_m = sigmac(epsilon,epsilon_c2,sigma_cd,n)*power(epsilon-epsilon_0,m);
end

function D_c = Dc(epsilon, epsilon_c2,sigma_cd,n)
if epsilon < 0
    D_c = 0;
else
    if (0 <= epsilon) && (epsilon <= epsilon_c2)
        D_c = sigma_cd*n*power(1-epsilon/epsilon_c2,n-1)/epsilon_c2;
    else
        D_c = 0;
    end
end
end

% Nao usada
% function D_s = Ds(epsilon_s,f_yd,epsilon_yd)
% if (epsilon_s < -epsilon_yd)
%     D_s = 0;
% else
%     if (-epsilon_yd < epsilon_s) && (epsilon_s < epsilon_yd)
%         D_s = f_yd/epsilon_yd*1000; %% ? pode estar errado
%     else
%         D_s = 0;
%     end
% end
% end

%% Esforco do concreto

function N_c = Nc(epsilon_0, k, epsilon_c2, sigma_cd,n,y_b,y_t,b,h,tol_k)
if abs(k*h) < tol_k
    N_c = sigmac(epsilon_0,epsilon_c2,sigma_cd,n)*b*h;
else
    epsilon_b = epsilon_0 + y_b*k;
    epsilon_t = epsilon_0 + y_t*k;
    N_c = b* DI0(epsilon_t,epsilon_b, epsilon_c2, sigma_cd,n)/k;
end
end


function M_c = Mc(epsilon_0, k, epsilon_c2, sigma_cd,n,y_b,y_t,b,h,tol_k)
if abs(k*h) < tol_k
    M_c = sigmac(epsilon_0,epsilon_c2,sigma_cd,n)*S(b,y_b,y_t);
else
    epsilon_b = epsilon_0 + y_b*k;
    epsilon_t = epsilon_0 + y_t*k;
    M_c = b* (DI1(epsilon_t,epsilon_b, epsilon_c2, sigma_cd,n)-epsilon_0*DI0(epsilon_t,epsilon_b, epsilon_c2, sigma_cd,n))/power(k,2);
end
end


function sigma_c = sigmac(epsilon,epsilon_c2,sigma_cd,n)
if epsilon < 0
    sigma_c = 0;
else
    if (0 <= epsilon) && (epsilon <= epsilon_c2)
        sigma_c = sigma_cd*(1-power(1-epsilon/epsilon_c2,n));
    else
        sigma_c = sigma_cd;
    end
end
end

%% Esforco do aco

function N_s = Ns(y_s,nb,phi,f_yd,epsilon_yd,epsilon_0,k)

A_si = pi*phi*phi*nb/4;
epsilon_s = epsilon_0 + y_s*k;
N_si = zeros(1,length(y_s));
for i=1:length(y_s)
    N_si(i) = A_si(i)*sigmas(epsilon_s(i),f_yd,epsilon_yd);
end

N_s = sum(N_si);
end

function M_s = Ms(y_s,nb,phi,f_yd,epsilon_yd,epsilon_0,k)

A_si = pi*phi*phi*nb/4;
epsilon_s = epsilon_0 + y_s*k;
M_si = zeros(1,length(y_s));
for i=1:length(y_s)
    M_si(i) = A_si(i)*sigmas(epsilon_s(i),f_yd,epsilon_yd)*y_s(i);
end

M_s = sum(M_si);
end

function sigma_s = sigmas(epsilon_s,f_yd,epsilon_yd)
if (epsilon_s < -epsilon_yd)
    sigma_s = -f_yd;
else
    if (-epsilon_yd < epsilon_s) && (epsilon_s < epsilon_yd)
        sigma_s = f_yd/epsilon_yd*(epsilon_s);
    else
        sigma_s = f_yd;
    end
end
end


%% Calculo da Matriz Jacobiana

function J_s = Js(E_s,epsilon_yd,phi,y_s,nb,epsilon_0,k)

A_si = pi*phi*phi*nb/4;
epsilon_s = epsilon_0 + y_s*k;

J_s = zeros(2);
for i=1:length(A_si)
if (epsilon_s(i) < -epsilon_yd)
    D_si = 0;
else
    if (-epsilon_yd <= epsilon_s(i)) && (epsilon_s(i) <= epsilon_yd)
        D_si = E_s; % MPa;
    else
        D_si = 0;
    end
end
J_s = J_s - D_si*A_si(i)*[1 , y_s(i); y_s(i) , power(y_s(i),2)];
end
end


function J_c = Jc(epsilon_t,epsilon_b, epsilon_0, epsilon_c2, sigma_cd,n,b,k,h,tol_k)
J_c = zeros(2);
if abs(k*h) < tol_k
    J_c(1,1) = -Dc(epsilon_0, epsilon_c2,sigma_cd,n)*b*h;
    J_c(1,2) = 0;
    J_c(2,1) = 0;
    J_c(2,2) = -Dc(epsilon_0, epsilon_c2,sigma_cd,n)*b*power(h,3)/12;
else
    J_c(1,1)=-b*DJm(epsilon_t,epsilon_b, epsilon_0, epsilon_c2, sigma_cd,n,0)/k;
    J_c(1,2)=-b*(DJm(epsilon_t,epsilon_b, epsilon_0, epsilon_c2, sigma_cd,n,1)-DI0(epsilon_t,epsilon_b, epsilon_c2, sigma_cd,n))/power(k,2);
    J_c(2,1)=-b*(DJm(epsilon_t,epsilon_b, epsilon_0, epsilon_c2, sigma_cd,n,1)-DI0(epsilon_t,epsilon_b, epsilon_c2, sigma_cd,n))/power(k,2);
    J_c(2,2)=-b*(DJm(epsilon_t,epsilon_b, epsilon_0, epsilon_c2, sigma_cd,n,2)-2*(DI1(epsilon_t,epsilon_b, epsilon_c2, sigma_cd,n)-epsilon_0*DI0(epsilon_t,epsilon_b, epsilon_c2, sigma_cd,n)))/power(k,3);
end
end

%% Verificacao do ELU

function ELU(epsilon_0,k,y_b,y_t,y_s,epsilon_c2,epsilon_cu)

epsilon_b = epsilon_0 + y_b*k;
epsilon_t = epsilon_0 + y_t*k;
epsilon_s = epsilon_0 + y_s*k;
epsilon_max = max(epsilon_b,epsilon_t);
epsilon_min = min(epsilon_b,epsilon_t);
epsilon_s_min = min(epsilon_s);

if epsilon_max > epsilon_cu
    fprintf("ELU ultrapassado: polo epsilon_cu!\n")
end
if epsilon_s_min < -10
    fprintf("ELU ultrapassado: polo epsilon_s!\n")
end
if epsilon_max - (epsilon_cu-epsilon_c2)*(epsilon_max-epsilon_min)/epsilon_cu > epsilon_c2
    fprintf("ELU ultrapassado: polo epsilon_c2!\n")
end
if (epsilon_max <= epsilon_cu) && (epsilon_s_min >= -10) && (epsilon_max - (epsilon_cu-epsilon_c2)*(epsilon_max-epsilon_min)/epsilon_cu <= epsilon_c2)
    fprintf("ELU ok!\n")
end
end

%% Calculo do momento estatico

function s = S(b,y_b,y_t)
s = b*(power(y_t,2)-power(y_b,2))/2; %secao retangular
end