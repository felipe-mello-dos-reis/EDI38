%% Problema classico de verificacao

% Limpando variaveis
clear all
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
epslon_yd = f_yd/E_s; % por mil
gamma_aco = 7850; %kg/m3


% Especificacoes de design
c = 0.02; % m
% diametros = [0.0125, 0.016, 0.025, 0.032];
% Diametros comerciais: 12,5 mm / 16 mm / 25 mm / 32 mm
phi_t = 0.025; % m
phi_c = 0.016; % m
N_t = 5; % Numero de barras de aco tracionadas
N_c = 3; % Numero de barras de aco comprimidas
b = 0.20; %m
h = 0.50; %m
A_st = N_t*pi*phi_t*phi_t/4;
A_sc = N_c*pi*phi_c*phi_c/4;
M_d = 316600*1e-6; %MN*m tinha que ser dado

d_t = c + phi_t/2;
d_c = c + phi_c/2;
%d_t = 0.05;
d = h-d_t;
epslon_c2 = 0; % por mil
epslon_cu = 0; % por mil
x_lim = 0;
n = 0;

tol = 1e-9;
i_max = 1e3;
% Especificacoes do problema
if fck <= 50
    x_lim = 0.45*d;
    n = 2;
    epslon_c2 = 2;
    epslon_cu = 3.5;
else
    if (50 < fck) && (fck < 90)
        x_lim = 0.35*d;
        n = 1.4 + 23.4*power((90-fck)/100,4);
        epslon_c2 = 2 + 0.085*power(fck-50,0.53);
        epslon_cu = 2.6 + 35*power((90-fck)/100,4);
    end
end

% Metodo da bisseccao
x_left = d_c;
x_right = d;
i = 0;
f_x = 1e6;
while ((x_right-x_left) > tol*d) && (i < i_max)  && (abs(f_x) > tol)
    i = i+1;
    f_x_left = F(x_left,sigma_cd,b,d,d_c,x_lim,n,epslon_c2,epslon_cu,epslon_yd,f_yd,A_st,A_sc);
    x = (x_right+x_left)/2; 
    f_x = F(x,sigma_cd,b,d,d_c,x_lim,n,epslon_c2,epslon_cu,epslon_yd,f_yd,A_st,A_sc);
    if (f_x*f_x_left < 0)
        x_right = x;
    else
        x_left = x;
    end
end
fprintf("x encontrado = %f\n", x)
fprintf("k_x encontrado = %f\n", x/d)
M_d_max = (Rcc(x,sigma_cd,b,d,x_lim,n,epslon_c2,epslon_cu)*d - Rcca(x,sigma_cd,b,d,x_lim,n,epslon_c2,epslon_cu) + Rsc(x,d,d_c,f_yd,epslon_c2,epslon_cu,epslon_yd,A_sc)*(d-d_c));
fprintf("M_d_max = %f\n", M_d_max) %sse ja eh o majorado, precisa dividir por 1.4
if (M_d_max - M_d > tol*M_d_max)
    fprintf("Verificacao OK!\n")
else
    if (abs(M_d-M_d_max) <= tol*M_d_max) &&(x < x_lim) % iminencia da ruptura
        fprintf("Verificacao OK!\n")
    else
        fprintf("Cessao rompe!\n")
    end
end

%% Calculo de R_cc

function R_cc = Rcc(x,sigma_cd,b,d,x_lim,n,epslon_c2,epslon_cu)
k_x = x/d;
x_lim2a2b = d*epslon_c2/(epslon_c2+10);
x_lim2b3 = d*epslon_cu/(epslon_cu+10);
if (0 <= x) && (x <= x_lim2a2b)
    R_cc = sigma_cd*b*d*(k_x-epslon_c2*(1-k_x)/(10*(n+1))*(1-power(1-10*k_x/(epslon_c2*(1-k_x)),n+1)));
else
    if (x_lim2a2b <= x) && (x <= x_lim2b3)
        R_cc = sigma_cd*b*d*(k_x-epslon_c2*(1-k_x)/(10*(n+1)));
    else
        if (x_lim2b3 < x) && (x <= d)
            R_cc = sigma_cd*b*d*k_x*(1-epslon_c2/(epslon_cu*(n+1)));
        else
            fprintf("x > d\n")
            R_cc = 0;
        end
    end
end
end

%% Calculo de R_cca

function R_cca = Rcca(x,sigma_cd,b,d,x_lim,n,epslon_c2,epslon_cu)
k_x = x/d;
x_lim2a2b = d*epslon_c2/(epslon_c2+10);
x_lim2b3 = d*epslon_cu/(epslon_cu+10);
if (0 <= x) && (x <= x_lim2a2b)
    sigma    R_cca = sigma_cd*b*d*d*(k_x*k_x/2-epslon_c2*(1-k_x)/(10*(n+1))*(k_x-epslon_c2*(1-k_x)/(10*(n+2))*(1-power(1-10*k_x/(epslon_c2*(1-k_x)),n+2))));
else
    if (x_lim2a2b <= x) && (x <= x_lim2b3)
        R_cca = sigma_cd*b*d*d*(k_x*k_x/2-epslon_c2*(1-k_x)/(10*(n+1))*(k_x-epslon_c2*(1-k_x)/(10*(n+2))));
    else
        if (x_lim2b3 < x) && (x <= d)
            R_cca = sigma_cd*b*d*d*k_x*k_x*(1/2-epslon_c2/(epslon_cu*(n+1))*(1-epslon_c2/(epslon_cu*(n+2))));
        else
            fprintf("x > d\n")
            R_cca = 0;
        end
    end
end
end

%% Funcao objetivo
function f = F(x,sigma_cd,b,d,d_c,x_lim,n,epslon_c2,epslon_cu,epslon_yd,f_yd,A_st,A_sc)
f = Rcc(x,sigma_cd,b,d,x_lim,n,epslon_c2,epslon_cu) + Rsc(x,d,d_c,f_yd,epslon_c2,epslon_cu,epslon_yd,A_sc) - Rst(x,d,f_yd,epslon_c2,epslon_cu,epslon_yd,A_st);
end

%% Calculo de R_st

function R_st = Rst(x,d,f_yd,epslon_c2,epslon_cu,epslon_yd,A_st)
x_lim34 = d*epslon_cu/(epslon_cu+epslon_yd);
if (0 <= x) && (x <= x_lim34)
    R_st = f_yd*A_st;
else
if (x_lim34 < x) && (x < d)
    R_st = A_st*f_yd/epslon_yd*(d-x)*epslon_cu/x;
else
    R_st = 0; % Erro -> preciso terminar ?
end
end
end

%% Calculo de R_sc

function R_sc = Rsc(x,d,d_c,f_yd,epslon_c2,epslon_cu,epslon_yd,A_sc)
x_ysc = (10*d_c+d*epslon_yd)/(epslon_yd+10);
x_lim2b3 = d*epslon_cu/(epslon_cu+10);
if (0 <= x) && (x <= x_lim2b3) 
    epslon_sc = (x-d_c)*10/(d-x);
else
if (x_lim2b3 < x) && (x <= d)
    epslon_sc = (x-d_c)*epslon_cu/x;
end
end
if epslon_sc <= epslon_yd
    sigma_sc = f_yd*epslon_sc/epslon_yd;
else
    sigma_sc = f_yd;
end
R_sc = sigma_sc*A_sc;
end

