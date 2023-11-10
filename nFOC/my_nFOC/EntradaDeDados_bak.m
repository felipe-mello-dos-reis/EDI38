function [Xc, Yc, INC, SIGMAcd, Ec2, Ecu, n, Xs, Ys, As, Nc, Ns, classe_aco, fyk, gamma_s, Es, fyd, Eyd, DEF, Nd, Mdx, Mdy, TOL_F, TOL_J, TOL_DEF, TOL_k, OBJETIVO_DA_ANALISE] = EntradaDeDados


%% COORDENADAS DOS PONTOS DA SEÇÃO POLIGONAL DE CONCRETO

% PONTO     Xc(cm)     Yc(cm)

Xc(1) = 0.0;
Yc(1) = 0.0;
Xc(2) = 20.0;
Yc(2) = 0.0;
Xc(3) = 20.0;
Yc(3) = 50.0;
Xc(4) = 0.0;
Yc(4) = 50.0;
Xc(5) = Xc(1);
Yc(5) = Yc(1);


%% INCIDENCIA DAS ARESTAS DA SEÇÃO POLIGONAL DE CONCRETO
% ARESTA        PONTO1          PONTO2
INC(1,1) = 1;
INC(2,1) = 2;
INC(1,2) = 2;
INC(2,2) = 3;
INC(1,3) = 3;
INC(2,3) = 4;
INC(1,4) = 4;
INC(2,4) = 1;
INC = INC';

%% PROPRIEDADES DO CONCRETO     
% Efeito Rüsch   fck (kN/cm2)        gamac    
Rusch = 0.85;
fck = 2.0;
gamma_c = 1.5;
SIGMAcd = Rusch*fck/gamma_c; %MPa

fck = fck*10; % convertendo para MPa

if fck <= 50
    % x_lim = 0.45*d;
    n = 2;
    Ec2 = 2;
    Ecu = 3.5;
else
    if (50 < fck) && (fck < 90)
        % x_lim = 0.35*d;
        n = 1.4 + 23.4*power((90-fck)/100,4);
        Ec2 = 2 + 0.085*power(fck-50,0.53);
        Ecu = 2.6 + 35*power((90-fck)/100,4);
    end
end

fck = fck/10; % retornando para kN/cm2

%% COORDENADAS DAS BARRAS DA ARMADURA E ÁREAS

Nc = 4;
Ns = [1, 1, 1, 1];
phi = 2.5;
% BARRA    Xs(cm)    Ys(cm)     Asi(cm²)

Xs(1) = 3.0;
Ys(1) = 3.0;
As(1) = Ns(1)*pi*phi^2/4;
Xs(2) = 17.0;
Ys(2) = 3.0;
As(2) = Ns(2)*pi*phi^2/4;
Xs(3) = 17.0;
Ys(3) = 47.0;
As(3) = Ns(3)*pi*phi^2/4;
Xs(4) = 3.0;
Ys(4) = 47.0;
As(4) = Ns(4)*pi*phi^2/4;

%% PROPRIEDADES DO AÇO     
% Classe   fYk(kN/cm²)     gamas     Es(kN/cm²)         
classe_aco = 'A';
fyk = 50.0;
gamma_s =  1.15;
Es = 20000.0;
fyd = fyk/gamma_s;
Eyd = fyd/Es*1000; % por mil



%%  DEFORMAÇÃO PRESCRITA PARA A SEÇÃO (CÁLCULO DOS ESFORÇOS RESISTÊNTES)     
% E0         kx(1/cm)      ky(1/cm)      
DEF(1) = 0.78947;
DEF(2) = -0.01486;
DEF(3) = -0.04710;

%% ESFORCOS DE CÁLCULO APLICADOS NA SEÇÃO (VERIFICAÇÃO DA SEÇÃO)    
% Nd(kN)            Mdx(kN.cm)      Mdy(kN.cm)      
Nd = 1000.0;
Mdx = -5000.0;
Mdy = -2000.0;  

% TOLERANCIA PARA A NORMA ADMENSIONAL DO VETOR F NO PROBLEMA DE VERIFICAÇÃO (MÉTODO DE NEWTON-RHAPSON)          
TOL_F = 1e-5 ;

% TOLERANCIA PARA O DETERMINANTE ADMENSIONAL DA MATRIZ J NO PROBLEMA DE VERIFICACAO (MÉTODO DE NEWTON-RHAPSON)                    
TOL_J = 1e-8;    

% PRECISAO PARA ÁREA DE ARMADURA NO PROBLEMA DE DIMENSIONAMENTO CONTINUO em cm2 (MÉTODO DA BISSEÇÃO)                     
%   1e-1 

% TOLERANCIAS ADMENSIONAIS PARA INTEDERMINACOES DO TIPO 1/DE OU 1/k
% DEFORMACOES(DE)       CURVATURAS(k.h)                    
TOL_DEF = 1e-3;
TOL_k = 1e-3;

% TOLERANCIA PARA A FLECHA ADMENSIONAL NO PROBLEMA DA VERIFACAÇÃO DA ESTABILIDADE DO PILAR (DIFERENÇAS FINITAS)                  
%     1e-5

% TOLERANCIA PARA O ESFORÇO NORMAL CRÍTICO ADMENSIONAL (MÉTODO DO PASSO)                  
%     1e-5

% NuMERO DE PONTOS POR TRECHO DOS DIAGRAMAS DE ESFORÇOS RESISTÊNTES              
%           500  

% NuMERO DE PONTOS POR CURVA MOMENTO-CURVATURA       
%           500  

% NuMERO DE PASSOS PARA A CONSTRUÇÃO DA TRAJETÓRIA DE EQUILÍBRIO DO PILAR              
%           500 

OBJETIVO_DA_ANALISE = 'ESFORCOS_RESISTENTES'; % Altere para a análise desejada

