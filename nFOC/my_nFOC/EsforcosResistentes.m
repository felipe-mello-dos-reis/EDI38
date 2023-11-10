function [ER, ERc, ERs] = EsforcosResistentes(Xc, Yc, INC, SIGMAcd, Ec2, Ecu, n, Xs, Ys, As, Nc, Ns, classe_aco, fyk, gamma_s, Es, fyd, Eyd, DEF, Nd, Mdx, Mdy, TOL_F, TOL_J, TOL_DEF, TOL_k, AREA, Sx, Sy, Ixx, Iyy, Ixy, SINAL_DA_CIRCUICAO, Ast, Tr_I, b, h, Ysmin, Ysmax)
% function [ER, SIGMAcd, Ec2, Ecu] = EsforcosResistentes(Xc, Yc, INC, SIGMAcd, Ec2, Ecu, n, Xs, Ys, As, Nc, Ns, classe_aco, fyk, gamma_s, Es, fyd, Eyd, DEF, Nd, Mdx, Mdy, TOL_F, TOL_J, TOL_DEF, TOL_k, AREA, Sx, Sy, Ixx, Iyy, Ixy, SINAL_DA_CIRCUICAO, Ast, Tr_I, b, h, Ysmin, Ysmax)


    % Parâmetros de entrada:
    % DEF - Vetor de deformações
    % h - Altura da seção
    % TOL_k - Tolerância para k
    % TOL_DEF - Tolerância para deformação
    % SINAL_DA_CIRCUICAO - Sinal de circuito
    % PROPRIEDADES_GEOMETRICAS_SECAO - Propriedades geométricas da seção
    % PROPRIEDADES_MATERIAIS - Propriedades dos materiais
    % GRANDEZAS_CINEMATICAS_SECAO - Grandezas cinemáticas da seção
    % ESFORCOS_RESISTENTES_E_DERIVADAS - Funções para cálculo de esforços resistentes e derivadas
    % ANALISE - Análise

    % Parâmetros definidos em PROPRIEDADES_GEOMETRICAS_SECAO, PROPRIEDADES_MATERIAIS, etc.

    % Inicialização das variáveis
    ERc = zeros(1, 3); % Esforços resistentes do concreto
    ERs = zeros(1, 3); % Esforços resistentes do aço


    % Grandezas cinemáticas
    DEF_1 = DEF(1);
    DEF_2 = DEF(2);
    DEF_3 = DEF(3);    

    % Chamar função para calcular tensão-deformação do concreto e obter SIGMA0
    ERc = EsforcosConcreto(DEF, h, TOL_k, TOL_DEF, Xc, Yc, INC, Xs, Ys, As, AREA, Sx, Sy, SIGMAcd, n, Ec2, Nc, SINAL_DA_CIRCUICAO);

    % Calcular esforços resistentes do aço
    for I = 1:Ns
        Ei = DEF(1) + DEF_3 * Xs(I) - DEF_2 * Ys(I);
        SIGMAi = TensaoDeformacaoAco(Ei, classe_aco, fyd, Eyd);
        ERs(1) = ERs(1) + SIGMAi * As(I);
        ERs(2) = ERs(2) - SIGMAi * Ys(I) * As(I);
        ERs(3) = ERs(3) + SIGMAi * Xs(I) * As(I);
    end

    % Esforços resistentes totais da seção
    ER = ERc + ERs;
end

% Função para cálculo dos esforços resistentes do concreto
function ERc = EsforcosConcreto(DEF, h, TOL_k, TOL_DEF, Xc, Yc, INC, Xs, Ys, As, AREA, Sx, Sy, SIGMAcd, n, Ec2, Nc, SINAL_DA_CIRCUICAO)
    if abs(DEF(2)*h) >= TOL_k
        ERc = zeros(1, 3);
        for I = 1:Nc
            Dxi = Xc(INC(I, 2)) - Xc(INC(I, 1));
            Dyi = Yc(INC(I, 2)) - Yc(INC(I, 1));
            Ei = DEF(1) + DEF(3)*Xc(INC(I, 1)) - DEF(2)*Yc(INC(I, 1));
            Ei_1 = DEF(1) + DEF(3)*Xc(INC(I, 2)) - DEF(2)*Yc(INC(I, 2));
            DEi = Ei_1 - Ei;
            hi = Xc(INC(I, 1))*Ei_1 - Xc(INC(I, 2))*Ei;
            gi = Yc(INC(I, 1))*Ei_1 - Yc(INC(I, 2))*Ei;
            [I1i, I2i, I3i, J1i, J2i, K1i] = IntegraisDiagramaTensaoDeformacao(Ei, SIGMAcd, n, Ec2);
            [I1i_1, I2i_1, I3i_1, J1i_1, J2i_1, K1i_1] = IntegraisDiagramaTensaoDeformacao(Ei_1, SIGMAcd, n, Ec2);
            DI2i = I2i_1 - I2i;
            DI3i = I3i_1 - I3i;
            DK1i = K1i_1 - K1i;
    
            if abs(DEi) <= TOL_DEF
                f1i = I1i;
                f2i = I1i*(Yc(INC(I, 1)) + Yc(INC(I, 2)))/2.0;
                f3i = f2i + I2i/DEF(2);
                f4i = I1i*(Xc(INC(I, 1)) + Xc(INC(I, 2)))/2.0;
            else
                f1i = DI2i/DEi;
                f2i = (gi*DI2i + Dyi*DK1i)/DEi^2.0;
                f3i = f2i + (1.0/DEF(2))*(DI3i/DEi);
                f4i = (hi*DI2i + Dxi*DK1i)/DEi^2.0;
            end
    
            if SINAL_DA_CIRCUICAO > 0
                ERc(1) = ERc(1) + Dxi*f1i;
                ERc(2) = ERc(2) - Dxi*f3i;
                ERc(3) = ERc(3) + Dxi*f4i;
            else
                ERc(1) = ERc(1) - Dxi*f1i;
                ERc(2) = ERc(2) + Dxi*f3i;
                ERc(3) = ERc(3) - Dxi*f4i;
            end
        end
        ERc = ERc*(1/DEF(2));
    elseif abs(DEF(3)*h) >= TOL_k
        ERc = zeros(1, 3);
        for I = 1:Nc
            Dxi = Xc(INC(I, 2)) - Xc(INC(I, 1));
            Dyi = Yc(INC(I, 2)) - Yc(INC(I, 1));
            Ei = DEF(1) + DEF(3)*Xc(INC(I, 1)) - DEF(2)*Yc(INC(I, 1));
            Ei_1 = DEF(1) + DEF(3)*Xc(INC(I, 2)) - DEF(2)*Yc(INC(I, 2));
            DEi = Ei_1 - Ei;
            hi = Xc(INC(I, 1))*Ei_1 - Xc(INC(I, 2))*Ei;
            gi = Yc(INC(I, 1))*Ei_1 - Yc(INC(I, 2))*Ei;
            [I1i, I2i, I3i, J1i, J2i, K1i] = IntegraisDiagramaTensaoDeformacao(Ei, SIGMAcd, n, Ec2);
            [I1i_1, I2i_1, I3i_1, J1i_1, J2i_1, K1i_1] = IntegraisDiagramaTensaoDeformacao(Ei_1, SIGMAcd, n, Ec2);
            DI2i = I2i_1 - I2i;
            DI3i = I3i_1 - I3i;
            DK1i = K1i_1 - K1i;
    
            if abs(DEi) <= TOL_DEF
                f1i = I1i;
                f2i = I1i*(Yc(INC(I, 1)) + Yc(INC(I, 2)))/2.0;
                f4i = I1i*(Xc(INC(I, 1)) + Xc(INC(I, 2)))/2.0;
                f5i = f4i - I2i/DEF(3);
            else
                f1i = DI2i/DEi;
                f2i = (gi*DI2i + Dyi*DK1i)/DEi^2.0;
                f4i = (hi*DI2i + Dxi*DK1i)/DEi^2.0;
                f5i = f4i - (1.0/DEF(3))*(DI3i/DEi);
            end
    
            if SINAL_DA_CIRCUICAO > 0
                ERc(1) = ERc(1) + Dyi*f1i;
                ERc(2) = ERc(2) - Dyi*f2i;
                ERc(3) = ERc(3) + Dyi*f5i;
            else
                ERc(1) = ERc(1) - Dyi*f1i;
                ERc(2) = ERc(2) + Dyi*f2i;
                ERc(3) = ERc(3) - Dyi*f5i;
            end
        end
        ERc = ERc*(1/DEF(3));
    elseif abs(DEF(2)*h) <= TOL_k && abs(DEF(3)*h) <= TOL_k
        SIGMA0 = TensaoDeformacaoConcreto(DEF(1), SIGMAcd, Ec2, n);
        ERc = [SIGMA0*AREA, SIGMA0*Sx, SIGMA0*Sy];
    end
    
end

% Função para calcular tensão-deformação do aço
function SIGMA = TensaoDeformacaoAco(E, classe_aco, Fyd, Eyd)
    switch classe_aco
        case 'A'
            if abs(E) <= Eyd
                SIGMA = Fyd * (E / Eyd);
            else
                if E < 0
                    SIGMA = -Fyd;
                elseif E == 0
                    SIGMA = 0;
                else
                    SIGMA = Fyd;
                end
            end
        case 'B'
            if abs(E) <= 0.7 * Eyd
                SIGMA = Fyd * (E / Eyd);
            elseif abs(E) > (Eyd + 2)
                if E < 0
                    SIGMA = -Fyd;
                elseif E == 0
                    SIGMA = 0;
                else
                    SIGMA = Fyd;
                end
            else
                if E < 0
                    SIGMA = -Fyd * (280 - 9 * Eyd + 3 * sqrt(800 * abs(E) + Eyd * (9 * Eyd - 560))) / 400;
                elseif E == 0
                    SIGMA = 0;
                else
                    SIGMA = Fyd * (280 - 9 * Eyd + 3 * sqrt(800 * abs(E) + Eyd * (9 * Eyd - 560))) / 400;
                end
            end
    end
end
