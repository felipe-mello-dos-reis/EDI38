%% Problema classico de dimensionamento
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
epslon_yd = f_yd/E_s; % por mil
gamma_aco = 7850; %kg/m3


% Especificacoes de design
c = 0.02; %m
b = 0.20; %m
h = 0.3905; %m
M_d = b*h*gamma_conc*10*12^2/2*1e-6*1.4; %MN*m #Majorar esforços
M_d = 196.8*1e-3; %MN*m #Majorar esforços

% diametros = [0.0125, 0.016, 0.025, 0.032]; %m
A_s_necessaria = zeros(2,1e3);

% for index=1:length(diametros)
    % phi = diametros(index); % Diametros comerciais: 12,5 mm / 16 mm / 25 mm / 32 mm
    % fprintf("Diametro = %f", phi)


    d_t = 0.0325;
    d_c = 0.0265; % Precisa inicializar em funçao de phi ou nao?
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
    % Aqui precisa do M_r_max?
    % M_r_max = Rcc(x_lim,sigma_cd,b,d,x_lim,n,epslon_c2,epslon_cu)*d - Rcca(x_lim,sigma_cd,b,d,x_lim,n,epslon_c2,epslon_cu) + sigmasc(x,d,d_c,f_yd,epslon_c2,epslon_cu,epslon_yd)*Asc(x,sigma_cd,b,d,d_c,f_yd,x_lim,n,epslon_c2,epslon_cu,epslon_yd,M_d)*(d-d_c);
    % if (M_d <= M_r_max)
        x_left = d_c;
        x_right = x_lim;
        i=0;
        f_x = 1e6;
        while ((x_right-x_left) > tol*d) && i < i_max
            i = i+1;
            x_m = (x_left + x_right)/2;
            f_x_m = F(x_m,sigma_cd,b,d,d_c,x_lim,n,f_yd,epslon_c2,epslon_cu,epslon_yd,M_d);
            A_s_necessaria(1,i) = f_x_m;
            A_s_necessaria(2,i) = x_m;
            x_1 = (x_left + x_m)/2;
            x_2 = (x_m + x_right)/2;
            f_x_1 = F(x_1,sigma_cd,b,d,d_c,x_lim,n,f_yd,epslon_c2,epslon_cu,epslon_yd,M_d);
            f_x_2 = F(x_2,sigma_cd,b,d,d_c,x_lim,n,f_yd,epslon_c2,epslon_cu,epslon_yd,M_d);
            if (f_x_1 > f_x_2)
                x_left = x_1;
            else
                x_right = x_2;
            end
        end
        fprintf("\tx_m encontrado = %f\n", x_m)
        fprintf("\tk_x_m encontrado = %f\n", x_m/d)
        fprintf("\tA_s encontrada (cm^2) = %f\n", A_s_necessaria(1,i)*1e4)
        % fprintf("\tNumero de barras = %f\n", A_st/(pi*phi*phi/4))
        % A_st_necessaria(index) = ceil(A_st/(pi*phi*phi/4))*pi*phi*phi/4;
        fprintf("\tA_st necessaria (cm^2) = %f\n", Ast(x_m,sigma_cd,b,d,d_c,f_yd,x_lim,n,epslon_c2,epslon_cu,epslon_yd,M_d)*1e4);
        fprintf("\tA_sc necessaria (cm^2) = %f\n", Asc(x_m,sigma_cd,b,d,d_c,f_yd,x_lim,n,epslon_c2,epslon_cu,epslon_yd,M_d)*1e4);

        % fprintf("\tNumero de barras necessarias = %d\n", ceil(A_st/(pi*phi*phi/4)))
        fprintf("\n")
    % else
        % fprintf("\tSecao rompe\n")
    % end

% end

%% Plotando a Area de aco em funcao de x
plot(A_s_necessaria(2,1:i)/d,A_s_necessaria(1,1:i)*1e4,"--r")
xlabel("k_x")
ylabel("A_{st}+A_{sc} (cm^2)")


%% Calculo de R_cc

function R_cc = Rcc(x,sigma_cd,b,d,x_lim,n,epslon_c2,epslon_cu) % output -> MN
k_x = x/d;
x_lim2a2b = d*epslon_c2/(epslon_c2+10);
x_lim2b3 = d*epslon_cu/(epslon_cu+10);
if (0 <= x) && (x <= x_lim2a2b)
    R_cc = sigma_cd*b*d*(k_x-epslon_c2*(1-k_x)/(10*(n+1))*(1-power(1-10*k_x/(epslon_c2*(1-k_x)),n+1)));
else
    if (x_lim2a2b <= x) && (x <= x_lim2b3)
        R_cc = sigma_cd*b*d*(k_x-epslon_c2*(1-k_x)/(10*(n+1)));
    else
        if (x_lim2b3 < x) && (x <= x_lim)
            R_cc = sigma_cd*b*d*k_x*(1-epslon_c2/(epslon_cu*(n+1)));
        else
            fprintf("\tx > x_lim\n")
            R_cc = 0;
        end
    end
end
end

%% Calculo de R_cca

function R_cca = Rcca(x,sigma_cd,b,d,x_lim,n,epslon_c2,epslon_cu) % output -> MNm
k_x = x/d;
x_lim2a2b = d*epslon_c2/(epslon_c2+10);
x_lim2b3 = d*epslon_cu/(epslon_cu+10);
if (0 <= x) && (x <= x_lim2a2b)
    R_cca = sigma_cd*b*d*d*(k_x*k_x/2-epslon_c2*(1-k_x)/(10*(n+1))*(k_x-epslon_c2*(1-k_x)/(10*(n+2))*(1-power(1-10*k_x/(epslon_c2*(1-k_x)),n+2))));
else
    if (x_lim2a2b <= x) && (x <= x_lim2b3)
        R_cca = sigma_cd*b*d*d*(k_x*k_x/2-epslon_c2*(1-k_x)/(10*(n+1))*(k_x-epslon_c2*(1-k_x)/(10*(n+2))));
    else
        if (x_lim2b3 < x) && (x <= x_lim)
            R_cca = sigma_cd*b*d*d*k_x*k_x*(1/2-epslon_c2/(epslon_cu*(n+1))*(1-epslon_c2/(epslon_cu*(n+2))));
        else
            fprintf("\tx > x_lim\n")
            R_cca = 0;
        end
    end
end
end

%% Funcao objetivo

function f = F(x,sigma_cd,b,d,d_c,x_lim,n,f_yd,epslon_c2,epslon_cu,epslon_yd,M_d) % output -> MNm
f = Ast(x,sigma_cd,b,d,d_c,f_yd,x_lim,n,epslon_c2,epslon_cu,epslon_yd,M_d) + Asc(x,sigma_cd,b,d,d_c,f_yd,x_lim,n,epslon_c2,epslon_cu,epslon_yd,M_d);
end

%% Area de aco comprimido

function A_sc = Asc(x,sigma_cd,b,d,d_c,f_yd,x_lim,n,epslon_c2,epslon_cu,epslon_yd,M_d)
A_sc = (M_d-Rcc(x,sigma_cd,b,d,x_lim,n,epslon_c2,epslon_cu)*d + Rcca(x,sigma_cd,b,d,x_lim,n,epslon_c2,epslon_cu))/(sigmasc(x,d,d_c,f_yd,epslon_c2,epslon_cu,epslon_yd)*(d-d_c));
if A_sc < 0
    fprintf("Erro!!! A_sc < 0\n")
end
end

%% Area de aco tracionado

function A_st = Ast(x,sigma_cd,b,d,d_c,f_yd,x_lim,n,epslon_c2,epslon_cu,epslon_yd,M_d)
A_st = Rcc(x,sigma_cd,b,d,x_lim,n,epslon_c2,epslon_cu)/sigmast(x,d,f_yd,epslon_c2,epslon_cu,epslon_yd) + (M_d-Rcc(x,sigma_cd,b,d,x_lim,n,epslon_c2,epslon_cu)*d + Rcca(x,sigma_cd,b,d,x_lim,n,epslon_c2,epslon_cu))/(sigmast(x,d,f_yd,epslon_c2,epslon_cu,epslon_yd)*(d-d_c));
if A_st < 0
    fprintf("Erro!!! A_st < 0\n")
end
end


function sigma_st = sigmast(x,d,f_yd,epslon_c2,epslon_cu,epslon_yd)
x_lim34 = d*epslon_cu/(epslon_cu+epslon_yd);
if (0 <= x) && (x <= x_lim34)
    sigma_st = f_yd;
else
    if (x_lim34 < x) && (x < d)
        sigma_st = f_yd/epslon_yd*(d-x)*epslon_cu/x;
    else
        sigma_st = 0; % Erro -> preciso terminar
    end
end
end

%% Calculo de sigma_sc

function sigma_sc = sigmasc(x,d,d_c,f_yd,epslon_c2,epslon_cu,epslon_yd)
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
end

