%% Problema classico de dimensionamento
% Limpando variaveis
clear all
clc

% Especificacoes concreto
fck = 40; %MPa
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
h = 0.40; %m
M_d = 1610*1e-6*1.4; %MN*m #Majorar esforços
diametros = [0.0125, 0.016, 0.020, 0.025, 0.032]; %m
A_st_necessaria = zeros(1,length(diametros));

for index=1:length(diametros)
    phi = diametros(index); % Diametros comerciais: 12,5 mm / 16 mm / 25 mm / 32 mm
    fprintf("Diametro = %f", phi)


    d_t = c + phi/2;
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
    M_r_max = Rcc(x_lim,sigma_cd,b,d,x_lim,n,epslon_c2,epslon_cu)*d - Rcca(x_lim,sigma_cd,b,d,x_lim,n,epslon_c2,epslon_cu);
    if (M_d <= M_r_max)
        x_left = 0;
        x_right = x_lim;
        i=0;
        f_x = 1e6;
        while ((x_right-x_left) > tol*d) && i < i_max && (abs(f_x) > 1e-9)
            i = i+1;
            f_x_left = F(x_left,sigma_cd,b,d,x_lim,n,epslon_c2,epslon_cu,M_d);
            x = (x_right+x_left)/2;
            f_x = F(x,sigma_cd,b,d,x_lim,n,epslon_c2,epslon_cu,M_d);
            if (f_x*f_x_left < 0)
                x_right = x;
            else
                x_left = x;
            end
        end
        fprintf("\tx encontrado = %f m\n", x)
        fprintf("\tk_x encontrado = %f\n", x/d)
        A_st = Rcc(x,sigma_cd,b,d,x_lim,n,epslon_c2,epslon_cu)/sigmast(x,d,f_yd,epslon_c2,epslon_cu,epslon_yd);
        fprintf("\tA_st encontrada = %f cm2\n", A_st*1e4)
        fprintf("\tNumero de barras = %f\n", A_st/(pi*phi*phi/4))
        A_st_necessaria(index) = ceil(A_st/(pi*phi*phi/4))*pi*phi*phi/4;
        fprintf("\tA_st necessaria = %f cm2\n", A_st_necessaria(index)*1e4);
        fprintf("\tNumero de barras necessarias = %d\n", ceil(A_st/(pi*phi*phi/4)))
        fprintf("\n")
    else
        fprintf("\tSecao rompe\n")
    end

end

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
function f = F(x,sigma_cd,b,d,x_lim,n,epslon_c2,epslon_cu,M_d) % output -> MNm
f = M_d - Rcc(x,sigma_cd,b,d,x_lim,n,epslon_c2,epslon_cu)*d + Rcca(x,sigma_cd,b,d,x_lim,n,epslon_c2,epslon_cu);
end
%% Area de aco
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
