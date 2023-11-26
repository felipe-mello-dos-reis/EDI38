function FlexaoCompostaObliqua

    TI = cputime;
    
    % Chama a função de entrada de dados
    [Xc, Yc, INC, SIGMAcd, Ec2, Ecu, n, Xs, Ys, As, Nc, Ns, classe_aco, fyk, gamma_s, Es, fyd, Eyd, DEF, Nd, Mdx, Mdy, TOL_F, TOL_J, TOL_DEF, TOL_k, OBJETIVO_DA_ANALISE] = EntradaDeDados;    
    % Traduz o sistema de coordenadas
    [Xc, Yc, Xs, Ys, Ysmin, Ysmax] = TranslacaoDoSistemaDeCoordenadas(Xc, Yc, Xs, Ys);
    [AREA, Sx, Sy, Ixx, Iyy, Ixy, SINAL_DA_CIRCUICAO, Ast, Tr_I, b, h, Ysmin, Ysmax] = PropriedadesGeometricasDaSecao(Xc, Yc, Xs, Ys, As, Nc, Ns);

    % Realiza a análise selecionada
    switch OBJETIVO_DA_ANALISE
        case 'ESFORCOS_RESISTENTES'
            [ER, ERc, ERs] = EsforcosResistentes(Xc, Yc, INC, SIGMAcd, Ec2, Ecu, n, Xs, Ys, As, Nc, Ns, classe_aco, fyk, gamma_s, Es, fyd, Eyd, DEF, Nd, Mdx, Mdy, TOL_F, TOL_J, TOL_DEF, TOL_k, AREA, Sx, Sy, Ixx, Iyy, Ixy, SINAL_DA_CIRCUICAO, Ast, Tr_I, b, h, Ysmin, Ysmax);
        case 'DIAGRAMAS_ESFORCOS_RESISTENTES'
             DiagramasEsforcosResistentes(Xc, Yc, INC, SIGMAcd, Ec2, Ecu, n, Xs, Ys, As, Nc, Ns, classe_aco, fyk, gamma_s, Es, fyd, Eyd, DEF, Nd, Mdx, Mdy, TOL_F, TOL_J, TOL_DEF, TOL_k, AREA, Sx, Sy, Ixx, Iyy, Ixy, SINAL_DA_CIRCUICAO, Ast, Tr_I, b, h, Ysmin, Ysmax);
        case 'VERIFICACAO_SECAO'
            [RUINA_SECAO, DEF, NORMA_F] = VerificacaoSecao(Xc, Yc, INC, SIGMAcd, Ec2, Ecu, n, Xs, Ys, As, Nc, Ns, classe_aco, fyk, gamma_s, Es, fyd, Eyd, DEF, Nd, Mdx, Mdy, TOL_F, TOL_J, TOL_DEF, TOL_k, AREA, Sx, Sy, Ixx, Iyy, Ixy, SINAL_DA_CIRCUICAO, Ast, Tr_I, b, h, Ysmin, Ysmax);
        case 'DIMENSIONAMENTO_CONTINUO'
            [DIMENSIONAVEL, Ast, DEF, NORMA_F] = DimensionamentoContinuo(Xc, Yc, INC, SIGMAcd, Ec2, Ecu, n, Xs, Ys, As, Nc, Ns, classe_aco, fyk, gamma_s, Es, fyd, Eyd, DEF, Nd, Mdx, Mdy, TOL_F, TOL_J, TOL_DEF, TOL_k, AREA, Sx, Sy, Ixx, Iyy, Ixy, SINAL_DA_CIRCUICAO, Ast, Tr_I, b, h, Ysmin, Ysmax);
        case 'DIMENSIONAMENTO_PRATICO'
            [DIMENSIONAVEL, phi, Ast, DEF, NORMA_F] = DimensionamentoPratico(Xc, Yc, INC, SIGMAcd, Ec2, Ecu, n, Xs, Ys, As, Nc, Ns, classe_aco, fyk, gamma_s, Es, fyd, Eyd, DEF, Nd, Mdx, Mdy, TOL_F, TOL_J, TOL_DEF, TOL_k, AREA, Sx, Sy, Ixx, Iyy, Ixy, SINAL_DA_CIRCUICAO, Ast, Tr_I, b, h, Ysmin, Ysmax);
        % case 'VERIFICACAO_ESTABILIDADE_PILAR'
        %     VerificacaoEstabilidadePilar;
        % case 'MOMENTO_CURVATURA'
        %     MomentoCurvatura;
        % case 'TRAJETORIA_EQUILIBRIO_PILAR'
        %     TrajetoriaEquilibrioPilar;
    end
    
    for I = 1:10
        fprintf('\n');
    end
    
    figure(3)
    LN = 1/DEF(2)*(DEF(1) + DEF(3)*linspace(min(min(Xc, Yc))*1.25, max(max(Xc, Yc))*1.25,100));
    plot(Xc, Yc, '-b', Xs, Ys, '*r',  linspace(min(min(Xc, Yc))*1.25, max(max(Xc, Yc))*1.25,100), LN, '--g')
    xlabel('x (cm)')
    ylabel('y (cm)')
    title('Secao Transversal')
    axis('equal')
    xlim([min(min(Xc, Yc)) max(max(Xc, Yc))]*1.25)
    ylim([min(min(Xc, Yc)) max(max(Xc, Yc))]*1.25)
    print -depsc2 ../../images/nFOC_q5.eps

    fprintf('SOLUÇÃO CONCLUÍDA\n\n');
    TF = cputime;
    fprintf('TRABALHO COMPUTACIONAL, SEGUNDOS %e\n', TF - TI);

end
