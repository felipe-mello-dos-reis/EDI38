function [RUINA_SECAO, DEF, NORMA_F] = VerificacaoSecao(Xc, Yc, INC, SIGMAcd, Ec2, Ecu, n, Xs, Ys, As, Nc, Ns, classe_aco, fyk, gamma_s, Es, fyd, Eyd, DEF, Nd, Mdx, Mdy, TOL_F, TOL_J, TOL_DEF, TOL_k, AREA, Sx, Sy, Ixx, Iyy, Ixy, SINAL_DA_CIRCUICAO, Ast, Tr_I, b, h, Ysmin, Ysmax);

    % Inicialização
    RUINA_SECAO = false;

    % Parâmetros
    TOL_NORMA_F = 1e-6;
    TOL_DET_J = 1e-6;

    % Inicialização de variáveis
    DEF = zeros(size(DEF));
    ER = EsforcosResistentes(Xc, Yc, INC, SIGMAcd, Ec2, Ecu, n, Xs, Ys, As, Nc, Ns, classe_aco, fyk, gamma_s, Es, fyd, Eyd, DEF, Nd, Mdx, Mdy, TOL_F, TOL_J, TOL_DEF, TOL_k, AREA, Sx, Sy, Ixx, Iyy, Ixy, SINAL_DA_CIRCUICAO, Ast, Tr_I, b, h, Ysmin, Ysmax);
    NORMA_F = norm([(Nd - ER(1))/(SIGMAcd * AREA), (Mdx - ER(2))/(SIGMAcd * AREA * h), (Mdy - ER(3))/(SIGMAcd * AREA * h)]);
    Nitemax_VERIFICACAO = 10000;
    I = 0;

    % Loop do método de Newton-Raphson
    while NORMA_F > TOL_NORMA_F
        I = I + 1;
        if I > Nitemax_VERIFICACAO
            RUINA_SECAO = true;
            break;
        end

        % Cálculo das derivadas dos esforços resistentes
        J = DerivadasDosEsforcosResistentes(Xc, Yc, INC, SIGMAcd, Ec2, Ecu, n, Xs, Ys, As, Nc, Ns, classe_aco, fyk, gamma_s, Es, fyd, Eyd, DEF, Nd, Mdx, Mdy, TOL_F, TOL_J, TOL_DEF, TOL_k, AREA, Sx, Sy, Ixx, Iyy, Ixy, SINAL_DA_CIRCUICAO, Ast, Tr_I, b, h, Ysmin, Ysmax);

        % Cálculo do determinante da matriz Jacobiana
        % DET_J = J(1,1)*J(2,2)*J(3,3)/(SIGMAcd^3*AREA*Tr_I^2) + J(2,1)*J(1,3)*J(3,2)/(SIGMAcd^3*(AREA*h)^2*Tr_I) + J(3,1)*J(1,2)*J(2,3)/(SIGMAcd^3*(AREA*h)^2*Tr_I) - J(2,3)*J(3,2)*J(1,1)/(SIGMAcd^3*AREA*Tr_I^2) - J(3,1)*J(1,3)*J(2,2)/(SIGMAcd^3*(AREA*h)^2*Tr_I) - J(2,1)*J(1,2)*J(3,3)/(SIGMAcd^3*(AREA*h)^2*Tr_I);


        if abs(det(J)/(SIGMAcd^3*AREA*Tr_I^2)) <= TOL_DET_J
            disp('MATRIZ JACOBIANA SINGULAR: RUINA DA SECAO TRANSVERSAL');
            RUINA_SECAO = true;
            break;
        else
            % Cálculo da variação das deformações
            % DDEF(1) = (J(2,2)*J(3,3) - J(2,3)*J(3,2))*(Nd - ER(1)) + (J(1,3)*J(3,2) - J(1,2)*J(3,3))*(Mdx - ER(2)) + (J(1,2)*J(2,3) - J(1,3)*J(2,2))*(Mdy - ER(3));
            % DDEF(2) = (J(2,3)*J(3,1) - J(2,1)*J(3,3))*(Nd - ER(1)) + (J(1,1)*J(3,3) - J(1,3)*J(3,1))*(Mdx - ER(2)) + (J(1,3)*J(2,1) - J(1,1)*J(2,3))*(Mdy - ER(3));
            % DDEF(3) = (J(2,1)*J(3,2) - J(2,2)*J(3,1))*(Nd - ER(1)) + (J(1,2)*J(3,1) - J(1,1)*J(3,2))*(Mdx - ER(2)) + (J(1,1)*J(2,2) - J(1,2)*J(2,1))*(Mdy - ER(3));
            
            % DET_J = J(1,1)*J(2,2)*J(3,3) + J(2,1)*J(1,3)*J(3,2) + J(3,1)*J(1,2)*J(2,3) - J(2,3)*J(3,2)*J(1,1) - J(3,1)*J(1,3)*J(2,2) - J(2,1)*J(1,2)*J(3,3);
            % DET_J = det(J);
            % DDEF = DDEF/DET_J;
            DDEF = J\(([Nd, Mdx, Mdy] - ER)');
            % Atualização das deformações
            DEF = DEF + DDEF';
        end

        % Recálculo dos esforços resistentes
        ER = EsforcosResistentes(Xc, Yc, INC, SIGMAcd, Ec2, Ecu, n, Xs, Ys, As, Nc, Ns, classe_aco, fyk, gamma_s, Es, fyd, Eyd, DEF, Nd, Mdx, Mdy, TOL_F, TOL_J, TOL_DEF, TOL_k, AREA, Sx, Sy, Ixx, Iyy, Ixy, SINAL_DA_CIRCUICAO, Ast, Tr_I, b, h, Ysmin, Ysmax);

        % Recálculo da norma da força
        NORMA_F = norm([(Nd - ER(1))/(SIGMAcd * AREA), (Mdx - ER(2))/(SIGMAcd * AREA * h), (Mdy - ER(3))/(SIGMAcd * AREA * h)]);

    end

    % Verificação ELU
    Emax = DEF(1) + DEF(3) * Xc(INC(1, 1)) - DEF(2) * Yc(INC(1, 1));
    Emin = DEF(1) + DEF(3) * Xc(INC(1, 1)) - DEF(2) * Yc(INC(1, 1));

    for I = 1:Nc
        Ei = DEF(1) + DEF(3) * Xc(INC(I, 1)) - DEF(2) * Yc(INC(I, 1));
        if Ei >= Emax
            Emax = Ei;
        end
        if Ei <= Emin
            Emin = Ei;
        end
    end

    Estar = Emax - ((Ecu - Ec2) / Ecu) * (Emax - Emin);

    if Ns == 0
        Esmin = Emin;
    else
        Esmin = DEF(1) + DEF(3) * Xs(1) - DEF(2) * Ys(1);

        for I = 1:Ns
            Ei = DEF(1) + DEF(3) * Xs(I) - DEF(2) * Ys(I);
            if Ei <= Esmin
                Esmin = Ei;
            end
        end
    end

    if Emax <= Ecu && Esmin >= -10 && Estar <= Ec2
        RUINA_SECAO = false;
    else
        RUINA_SECAO = true;
        disp('FALHA NA VERIFICACAO DO ELU');
    end
end
