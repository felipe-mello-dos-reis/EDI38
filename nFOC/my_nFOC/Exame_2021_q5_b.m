function [Xc, Yc, INC, SIGMAcd, Ec2, Ecu, n, Xs, Ys, As, Nc, Ns, classe_aco, fyk, gamma_s, Es, fyd, Eyd, DEF, Nd, Mdx, Mdy, TOL_F, TOL_J, TOL_DEF, TOL_k, OBJETIVO_DA_ANALISE] = EntradaDeDados

    
    %% COORDENADAS DOS PONTOS DA SEÇÃO POLIGONAL DE CONCRETO
    
    % PONTO     Xc(cm)     Yc(cm)
    
    % Define the data
    Xc = [-15,15,15,7.5,7.5,15,15,-15,-15,-7.5,-7.5,-15];
    Yc = [-40,-40,-30,-25,25,30,40,40,30,25,-25,-30];
    Xc(length(Xc) + 1) = Xc(1);
    Yc(length(Yc) + 1) = Yc(1);

    Nc = length(Xc)-1;
    INC = zeros(2, Nc);
    
    
    %% INCIDENCIA DAS ARESTAS DA SEÇÃO POLIGONAL DE CONCRETO
    % ARESTA        PONTO1          PONTO2
    % Preenche as arestas com base nas coordenadas dos vértices
    for i = 1:Nc
        INC(1, i) = i;
        if i == Nc
            INC(2, i) = 1;
        else
            INC(2, i) = i + 1;
        end
    end

    INC = INC';
    
    %% PROPRIEDADES DO CONCRETO     
    % Efeito Rüsch   fck (kN/cm2)        gamac    
    Rusch = 0.85;
    fck = 6.5;
    gamma_c = 1.4;
    SIGMAcd = Rusch*fck/gamma_c; %MPa
    
    fck = fck*10; % convertendo para MPa
    
    if fck <= 50
        % x_lim = 0.45*d;
        n = 2;
        Ec2 = 2;
        Ecu = 3.5;
    else
        if (50 < fck) && (fck <= 90)
            % x_lim = 0.35*d;
            n = 1.4 + 23.4*power((90-fck)/100,4);
            Ec2 = 2 + 0.085*power(fck-50,0.53);
            Ecu = 2.6 + 35*power((90-fck)/100,4);
        end
    end
    
    fck = fck/10; % retornando para kN/cm2
    
    %% COORDENADAS DAS BARRAS DA ARMADURA E ÁREAS
    % BARRA    Xs(cm)    Ys(cm)     Asi(cm²)
    
    % Define the data
    Xs = [-12,0,12,-12,12,-12,12,-12,0,12];
    Ys = [-37,-37,-37,-33,-33,33,33,37,37,37];

    phi = 2.5;
    As = pi()*phi^2/4*ones(size(Xs));
    Ns = length(Xs);
    
    %% PROPRIEDADES DO AÇO     
    % Classe   fYk(kN/cm²)     gamas     Es(kN/cm²)         
    classe_aco = 'A';
    fyk = 50.0;
    gamma_s =  1.15;
    Es = 21000.0;
    fyd = fyk/gamma_s;
    Eyd = fyd/Es*1000; % por mil
    
    
    
    %%  DEFORMAÇÃO PRESCRITA PARA A SEÇÃO (CÁLCULO DOS ESFORÇOS RESISTÊNTES)     
    % E0         kx(1/cm)      ky(1/cm)      
    DEF(1) = -0.30407777;
    DEF(2) = 0.030407777;
    DEF(3) = -0.12163111;
    DEF(2) = Ecu/90;
    DEF(1) = -10*DEF(2);
    DEF(3) = -4*DEF(2)
    
    %% ESFORCOS DE CÁLCULO APLICADOS NA SEÇÃO (VERIFICAÇÃO DA SEÇÃO)    
    % Nd(kN)            Mdx(kN.cm)      Mdy(kN.cm)      

    Nd = 2000;
    Mdx = -55000;
    Mdy = 30000;  

    % Nd = 1000.0;
    % Mdx = -14000.0;
    % Mdy = -2000.0;

    % Nd = 1000.0;
    % Mdx = -16244.0;
    % Mdy = -2000.0;

    % Nd = 1000.0;
    % Mdx = -5000.0;
    % Mdy = -2000.0;  

    % Nd = 1000.0;
    % Mdx = -3000.0;
    % Mdy = -8000.0;  
    
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
    
    OBJETIVO_DA_ANALISE = 'DIMENSIONAMENTO_CONTINUO'; % Altere para a análise desejada
    
    