%% Nova Flexao Normal Composta

% Limpando variaveis
clear all
close all
clc

% Especificacoes concreto
fck = 30; %MPa
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
c = 0.02; % m
b = 0.20; % m
h = 0.32; % m
d = 0.05;
nc = 3;
nb = [2 2 2];
phi = 0.010; % m
y_s = linspace(-(h/2-d),(h/2-d),nc)
% M_d = 110000*1e-6*1.4; %MN*m #Majorar esforços


epsilon_c2 = 0; % por mil
epsilon_cu = 0; % por mil
x_lim = 0;
n = 0;

tol = 1e-9;

% Especificacoes do problema
if fck <= 50
    x_lim = 0.45*d;
    n = 2;
    epsilon_c2 = 2;
    epsilon_cu = 3.5;
else
    if (50 < fck) && (fck < 90)
        x_lim = 0.35*d;
        n = 1.4 + 23.4*power((90-fck)/100,4);
        epsilon_c2 = 2 + 0.085*power(fck-50,0.53);
        epsilon_cu = 2.6 + 35*power((90-fck)/100,4);
    end
end

%% Plotando os esforços resistentes

y_t=h/2;
y_b=-h/2;
%plot_epsilon_0_k(epsilon_c2,epsilon_cu,h,y_t,y_b,y_s)
%plot_N_r_M_r(epsilon_c2,epsilon_cu,sigma_cd,n,b,f_yd,epsilon_yd,h,y_t,y_b,y_s,phi,nb)

%% Testando um valor especifico
epsilon_0 = 0.5;
k = 9.375;

A_si = pi*phi*phi*nb/4;
epsilon_s = epsilon_0 + y_s*k;
N_si = zeros(1,length(y_s));
M_si = zeros(1,length(y_s));
for i=1:length(y_s)
    N_si(i) = A_si(i)*sigmas(epsilon_s(i),f_yd,epsilon_yd);
    M_si(i) = A_si(i)*sigmas(epsilon_s(i),f_yd,epsilon_yd)*y_s(i);
end

N_s = sum(N_si);
M_s = sum(M_si);
N_c = Nc(epsilon_0, k, epsilon_c2, sigma_cd,n,y_b,y_t,b,h);
M_c = Mc(epsilon_0, k, epsilon_c2, sigma_cd,n,y_b,y_t,b);
N_r = double(N_c + N_s)
M_r = double(M_c + M_s)

function I_m = Im(epsilon, epsilon_c2, sigma_cd,n,m)
syms csi
if epsilon <= 0
    I_m = 0;
else
    if (0 < epsilon) && (epsilon <= epsilon_c2)
        I_m = int(sigma_cd*(1-power(1-csi/epsilon_c2,n))*power(csi,m),csi,0,epsilon);
    else
        if epsilon > epsilon_c2
            I_m = int(sigma_cd*(1-power(1-csi/epsilon_c2,n))*power(csi,m),csi,0,epsilon_c2) + sigma_cd*int(power(csi,m),csi,epsilon_c2,epsilon);
        end
    end
end
end

function plot_epsilon_0_k(epsilon_c2,epsilon_cu,h,y_t,y_b,y_s)
c = epsilon_c2;
u = epsilon_cu;
A = [0;c];
B = [u/h;u*y_t/h];
C = [(u+10)/(y_t-min(y_s));u-y_t*((u+10)/(y_t-min(y_s)))];
D = [0;-10];
E = [(u+10)/(y_b-max(y_s));u-y_b*((u+10)/(y_b-max(y_s)))];
F = [-u/h;-u*y_b/h];
points = [A B C D E F A]
figure(1);
set(gca, 'FontSize', 12);
set(gca, 'TickLabelInterpreter', 'latex');
plot(points(1,:),points(2,:),'-b', 'LineWidth', 2);
% Get the current axes
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';xlabel('$\kappa \; (1/m)$', Interpreter='latex')
ylabel('$\epsilon_0$', Interpreter='latex')
xlim(1.1*[min(points(1,:)) max(points(1,:))])
ylim(1.1*[min(points(2,:)) max(points(2,:))])
title('Regiao viavel para o dimensionamento', Interpreter='latex')

end

function plot_N_r_M_r(epsilon_c2,epsilon_cu,sigma_cd,n,b,f_yd,epsilon_yd,h,y_t,y_b,y_s,phi,nb)
c = epsilon_c2;
u = epsilon_cu;
A = [0;c];
B = [u/h;u*y_t/h];
C = [(u+10)/(y_t-min(y_s));u-y_t*((u+10)/(y_t-min(y_s)))];
D = [0;-10];
E = [(u+10)/(y_b-max(y_s));u-y_b*((u+10)/(y_b-max(y_s)))];
F = [-u/h;-u*y_b/h];
points = [A B C D E F A];
% Number of points you want to generate
num_points = 10;
kappa_values = [];
epsilon_0_values = [];
for i=1:6
    kappa_values = [kappa_values, linspace(points(1,i), points(1,i+1), num_points)];
    epsilon_0_values = [epsilon_0_values,  linspace(points(2,i), points(2,i+1), num_points)];

end
N_r_values = zeros(1,length(kappa_values));
M_r_values = zeros(1,length(epsilon_0_values));

for j=1:length(kappa_values)
    fprintf("Iteracao %d\n",j)
    N_c = Nc(epsilon_0_values(j), kappa_values(j), epsilon_c2, sigma_cd,n,y_b,y_t,b,h);
    M_c = Mc(epsilon_0_values(j), kappa_values(j), epsilon_c2, sigma_cd,n,y_b,y_t,b);
    A_si = pi*phi*phi*nb/4;
    epsilon_s = epsilon_0_values(j) + y_s*kappa_values(j);
    N_si = zeros(1,length(y_s));
    M_si = zeros(1,length(y_s));
    for i=1:length(y_s)
        N_si(i) = A_si(i)*sigmas(epsilon_s(i),f_yd,epsilon_yd);
        M_si(i) = A_si(i)*sigmas(epsilon_s(i),f_yd,epsilon_yd)*y_s(i);
    end
    N_s = sum(N_si);
    M_s = sum(M_si);
    N_r_values(j) = (N_c + N_s)*1e3;
    M_r_values(j) = (M_c + M_s)*1e3;
end

figure(2);
plot(N_r_values, M_r_values, '-r', 'LineWidth', 2);
% Get the current axes
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

function DI_m = DIm(epsilon_t,epsilon_b, epsilon_c2, sigma_cd,n,m)
DI_m = Im(epsilon_t, epsilon_c2, sigma_cd,n,m)-Im(epsilon_b, epsilon_c2, sigma_cd,n,m);
end

function N_c = Nc(epsilon_0, k, epsilon_c2, sigma_cd,n,y_b,y_t,b,h)
tol = 1e-10;
if abs(k) < tol
    N_c = sigmac(epsilon_0,epsilon_c2,sigma_cd,n)*b*h;
else
    epsilon_b = epsilon_0 + y_b*k;
    epsilon_t = epsilon_0 + y_t*k;
    N_c = b* DIm(epsilon_t,epsilon_b, epsilon_c2, sigma_cd,n,0)/k;
end
end


function M_c = Mc(epsilon_0, k, epsilon_c2, sigma_cd,n,y_b,y_t,b)
tol = 1e-10;
if abs(k) < tol
    M_c = sigmac(epsilon_0,epsilon_c2,sigma_cd,n)*S(b,y_b,y_t);
else
    epsilon_b = epsilon_0 + y_b*k;
    epsilon_t = epsilon_0 + y_t*k;
    M_c = b* (DIm(epsilon_t,epsilon_b, epsilon_c2, sigma_cd,n,1)-epsilon_0*DIm(epsilon_t,epsilon_b, epsilon_c2, sigma_cd,n,0))/power(k,2);
end
end


function sigma_c = sigmac(epsilon,epsilon_c2,sigma_cd,n)
if epsilon <= 0
    sigma_c = 0;
else
    if (0 < epsilon) && (epsilon <= epsilon_c2)
        sigma_c = sigma_cd*(1-power(1-epsilon/epsilon_c2,n));
    else
        sigma_c = sigma_cd;
    end
end
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

function s = S(b,y_b,y_t)
syms y
s = double(b*int(y,y,y_b,y_t));
end