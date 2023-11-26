function DiagramasEsforcosResistentes(Xc, Yc, INC, SIGMAcd, Ec2, Ecu, n, Xs, Ys, As, Nc, Ns, classe_aco, fyk, gamma_s, Es, fyd, Eyd, DEF, Nd, Mdx, Mdy, TOL_F, TOL_J, TOL_DEF, TOL_k, AREA, Sx, Sy, Ixx, Iyy, Ixy, SINAL_DA_CIRCUICAO, Ast, Tr_I, b, h, Ysmin, Ysmax)
    % Define os pontos PA, PB, PC, PD, PE e PF
    Ymax = max(Yc);
    Ymin = min(Yc);
    Ysmax = max(Ys);
    Ysmin = min(Ys);
    PA = [0; Ec2];
    PB = [Ecu / h; Ecu * Ymax / h];
    PC = [(Ecu + 10) / (Ymax - Ysmin); Ecu - Ymax * (Ecu + 10) / (Ymax - Ysmin)];
    PD = [0; -10];
    PE = [(Ecu + 10) / (Ymin - Ysmax); Ecu - Ymin * (Ecu + 10) / (Ymin - Ysmax)];
    PF = [-Ecu / h; -Ecu * Ymin / h];

    % Inicializa a matriz de resultados dos diagramas
    NPONTOS_DIAGRAMA = 100; % Defina o número de pontos desejado
    DIAGRAMA_Nr_Mrx = zeros(6 * NPONTOS_DIAGRAMA, 2);

    % Define a deflexão como zero
    % DEF(3) = 0.0;

    % Inicializa o índice para preencher os resultados
    II = 0;

    % Loop para preencher os diagramas de Nr e Mrx para as arestas AB, BC, CD, DE, EF e FA
    for I = 1:NPONTOS_DIAGRAMA
        % Calcula o ponto intermediário P para a aresta AB
        P = PA + (PB - PA) * I / NPONTOS_DIAGRAMA;
        DEF(1) = P(2);
        DEF(2) = P(1);
        [ER, ~] = EsforcosResistentes(Xc, Yc, INC, SIGMAcd, Ec2, Ecu, n, Xs, Ys, As, Nc, Ns, classe_aco, fyk, gamma_s, Es, fyd, Eyd, DEF, Nd, Mdx, Mdy, TOL_F, TOL_J, TOL_DEF, TOL_k, AREA, Sx, Sy, Ixx, Iyy, Ixy, SINAL_DA_CIRCUICAO, Ast, Tr_I, b, h, Ysmin, Ysmax); % Substitua com a função correta
        II = II + 1;
        DIAGRAMA_Nr_Mrx(II, 1) = ER(1);
        DIAGRAMA_Nr_Mrx(II, 2) = ER(2);
    end

    % Repita o loop acima para as outras arestas (BC, CD, DE, EF, FA)

    for I = 1:NPONTOS_DIAGRAMA
        % Calcula o ponto intermediário P para a aresta BC
        P = PB + (PC - PB) * I / NPONTOS_DIAGRAMA;
        DEF(1) = P(2);
        DEF(2) = P(1);
        [ER, ~] = EsforcosResistentes(Xc, Yc, INC, SIGMAcd, Ec2, Ecu, n, Xs, Ys, As, Nc, Ns, classe_aco, fyk, gamma_s, Es, fyd, Eyd, DEF, Nd, Mdx, Mdy, TOL_F, TOL_J, TOL_DEF, TOL_k, AREA, Sx, Sy, Ixx, Iyy, Ixy, SINAL_DA_CIRCUICAO, Ast, Tr_I, b, h, Ysmin, Ysmax);
        II = II + 1;
        DIAGRAMA_Nr_Mrx(II, 1) = ER(1);
        DIAGRAMA_Nr_Mrx(II, 2) = ER(2);
    end

    % Exemplo para a aresta CD
    for I = 1:NPONTOS_DIAGRAMA
        % Calcula o ponto intermediário P para a aresta CD
        P = PC + (PD - PC) * I / NPONTOS_DIAGRAMA;
        DEF(1) = P(2);
        DEF(2) = P(1);
        [ER, ~] = EsforcosResistentes(Xc, Yc, INC, SIGMAcd, Ec2, Ecu, n, Xs, Ys, As, Nc, Ns, classe_aco, fyk, gamma_s, Es, fyd, Eyd, DEF, Nd, Mdx, Mdy, TOL_F, TOL_J, TOL_DEF, TOL_k, AREA, Sx, Sy, Ixx, Iyy, Ixy, SINAL_DA_CIRCUICAO, Ast, Tr_I, b, h, Ysmin, Ysmax);
        II = II + 1;
        DIAGRAMA_Nr_Mrx(II, 1) = ER(1);
        DIAGRAMA_Nr_Mrx(II, 2) = ER(2);
    end

    % Repita o loop para a aresta DE
    for I = 1:NPONTOS_DIAGRAMA
        % Calcula o ponto intermediário P para a aresta DE
        P = PD + (PE - PD) * I / NPONTOS_DIAGRAMA;
        DEF(1) = P(2);
        DEF(2) = P(1);
        [ER, ~] = EsforcosResistentes(Xc, Yc, INC, SIGMAcd, Ec2, Ecu, n, Xs, Ys, As, Nc, Ns, classe_aco, fyk, gamma_s, Es, fyd, Eyd, DEF, Nd, Mdx, Mdy, TOL_F, TOL_J, TOL_DEF, TOL_k, AREA, Sx, Sy, Ixx, Iyy, Ixy, SINAL_DA_CIRCUICAO, Ast, Tr_I, b, h, Ysmin, Ysmax);
        II = II + 1;
        DIAGRAMA_Nr_Mrx(II, 1) = ER(1);
        DIAGRAMA_Nr_Mrx(II, 2) = ER(2);
    end

    % Repita o loop para a aresta EF
    for I = 1:NPONTOS_DIAGRAMA
        % Calcula o ponto intermediário P para a aresta EF
        P = PE + (PF - PE) * I / NPONTOS_DIAGRAMA;
        DEF(1) = P(2);
        DEF(2) = P(1);
        [ER, ~] = EsforcosResistentes(Xc, Yc, INC, SIGMAcd, Ec2, Ecu, n, Xs, Ys, As, Nc, Ns, classe_aco, fyk, gamma_s, Es, fyd, Eyd, DEF, Nd, Mdx, Mdy, TOL_F, TOL_J, TOL_DEF, TOL_k, AREA, Sx, Sy, Ixx, Iyy, Ixy, SINAL_DA_CIRCUICAO, Ast, Tr_I, b, h, Ysmin, Ysmax);
        II = II + 1;
        DIAGRAMA_Nr_Mrx(II, 1) = ER(1);
        DIAGRAMA_Nr_Mrx(II, 2) = ER(2);
    end

    % Repita o loop para a aresta FA
    for I = 1:NPONTOS_DIAGRAMA
        % Calcula o ponto intermediário P para a aresta FA
        P = PF + (PA - PF) * I / NPONTOS_DIAGRAMA;
        DEF(1) = P(2);
        DEF(2) = P(1);
        [ER, ~] = EsforcosResistentes(Xc, Yc, INC, SIGMAcd, Ec2, Ecu, n, Xs, Ys, As, Nc, Ns, classe_aco, fyk, gamma_s, Es, fyd, Eyd, DEF, Nd, Mdx, Mdy, TOL_F, TOL_J, TOL_DEF, TOL_k, AREA, Sx, Sy, Ixx, Iyy, Ixy, SINAL_DA_CIRCUICAO, Ast, Tr_I, b, h, Ysmin, Ysmax);
        II = II + 1;
        DIAGRAMA_Nr_Mrx(II, 1) = ER(1);
        DIAGRAMA_Nr_Mrx(II, 2) = ER(2);
    end





    % Crie um arquivo de saída e escreva os resultados
    fid = fopen('DIAGRAMAS_ESFORCOS_RESISTENTES.txt', 'w');
    fprintf(fid, '******************** DIAGRAMAS DE ESFORCOS RESISTENTES ********************\n');
    fprintf(fid, 'PONTOS DO DIAGRAMA Nr Mrx\n\n');
    fprintf(fid, '   PONTO          Nr               Mrx\n');
    
    for I = 1:6 * NPONTOS_DIAGRAMA
        fprintf(fid, '%5d     %13.6E     %13.6E\n', I, DIAGRAMA_Nr_Mrx(I, 1), DIAGRAMA_Nr_Mrx(I, 2));
    end
    points = [PA PB PC PD PE PF PA];
    figure(1);
    set(gca, 'FontSize', 12);
    plot(points(1,:),points(2,:),'-b', 'LineWidth', 2);
    xlabel('\kappa_x (1/m)')
    ylabel('\epsilon_0')
    xlim(1.1*[min(points(1,:)) max(points(1,:))])
    ylim(1.1*[min(points(2,:)) max(points(2,:))])
    title('Regiao viavel para o dimensionamento')
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    print -depsc2 ../../images/nFOC_q4_c_1.eps


    figure(2);
    set(gca, 'FontSize', 12);    
    plot(DIAGRAMA_Nr_Mrx(:, 1), DIAGRAMA_Nr_Mrx(:, 2), '-r', 'LineWidth', 2);
    xlabel('N_r (kN)');
    ylabel('M_{rx} (kNcm)');
    xlim([min(DIAGRAMA_Nr_Mrx(:, 1))*1.1, max(DIAGRAMA_Nr_Mrx(:, 1))*1.1]);
    ylim([min(DIAGRAMA_Nr_Mrx(:, 2))*1.1, max(DIAGRAMA_Nr_Mrx(:, 2))*1.1]);
    title('Esforcos Resistentes (\kappa_y = 0)');
    ax = gca;
    ax.XAxisLocation = 'origin';
    ax.YAxisLocation = 'origin';
    print -depsc2 ../../images/nFOC_q4_c_2.eps

    fclose(fid);
end 
