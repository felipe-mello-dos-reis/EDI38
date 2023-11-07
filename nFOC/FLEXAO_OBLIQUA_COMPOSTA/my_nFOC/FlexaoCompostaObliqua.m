function FlexaoCompostaObliqua

    TI = cputime;
    
    % Chama a função de entrada de dados
    [Xc, Yc, INC, SIGMAcd, Ec2, Ecu, n, Xs, Ys, As, Nc, Ns, classe_aco, fyk, gamma_s, Es, fyd, Eyd, DEF, Nd, Mdx, Mdy, TOL_F, TOL_J, TOL_DEF, TOL_k, OBJETIVO_DA_ANALISE] = EntradaDeDados;    
    % Traduz o sistema de coordenadas
    [AREA, Sx, Sy, Ixx, Iyy, Ixy, SINAL_DA_CIRCUICAO, Ast, Tr_I, b, h, Ysmin, Ysmax] = PropriedadesGeometricasDaSecao(Xc, Yc, Xs, Ys, As, Nc, Ns);

    % Realiza a análise selecionada
    switch OBJETIVO_DA_ANALISE
        case 'ESFORCOS_RESISTENTES'
            [ER, SIGMAcd, Ec2, Ecu] = EsforcosResistentes(Xc, Yc, INC, SIGMAcd, Ec2, Ecu, n, Xs, Ys, As, Nc, Ns, classe_aco, fyk, gamma_s, Es, fyd, Eyd, DEF, Nd, Mdx, Mdy, TOL_F, TOL_J, TOL_DEF, TOL_k, AREA, Sx, Sy, Ixx, Iyy, Ixy, SINAL_DA_CIRCUICAO, Ast, Tr_I, b, h, Ysmin, Ysmax);
        case 'DIAGRAMAS_ESFORCOS_RESISTENTES'
            [Mdxs, Mdys, Nds] = DiagramasEsforcosResistentes(Xc, Yc, Xs, Ys, As, Nc, Ns);
        case 'VERIFICACAO_SECAO'
            [Nd, Mdx, Mdy] = VerificacaoSecao(Nds, Mdxs, Mdys);
        case 'DIMENSIONAMENTO_CONTINUO'
            DimensionamentoContinuo;
        case 'DIMENSIONAMENTO_PRATICO'
            DimensionamentoPratico;
        case 'VERIFICACAO_ESTABILIDADE_PILAR'
            VerificacaoEstabilidadePilar;
        case 'MOMENTO_CURVATURA'
            MomentoCurvatura;
        case 'TRAJETORIA_EQUILIBRIO_PILAR'
            TrajetoriaEquilibrioPilar;
    end
    
    for I = 1:5
        fprintf('\n');
    end
    
    fprintf('SOLUÇÃO CONCLUÍDA\n\n');
    TF = cputime;
    fprintf('TRABALHO COMPUTACIONAL, SEGUNDOS %e\n', TF - TI);
end
