function [DIMENSIONAVEL, Ast, DEF, NORMA_F] = DimensionamentoContinuo(Xc, Yc, INC, SIGMAcd, Ec2, Ecu, n, Xs, Ys, As, Nc, Ns, classe_aco, fyk, gamma_s, Es, fyd, Eyd, DEF, Nd, Mdx, Mdy, TOL_F, TOL_J, TOL_DEF, TOL_k, AREA, Sx, Sy, Ixx, Iyy, Ixy, SINAL_DA_CIRCUICAO, Ast, Tr_I, b, h, Ysmin, Ysmax)

    DIMENSIONAVEL = false;

    As_left = 0.0;
    As_right = 0.04*AREA;
    It_max = 10000;
    I = 0;
    psi = 1/2;
    TOL_As = 1e-6;

    % Iteracao zero

    Ast = As_left;
    phi = sqrt(Ast*4/(pi()*Ns));
    As = pi()*phi^2/4*ones(size(Xs));
    [RUINA_SECAO, DEF, NORMA_F] = VerificacaoSecao(Xc, Yc, INC, SIGMAcd, Ec2, Ecu, n, Xs, Ys, As, Nc, Ns, classe_aco, fyk, gamma_s, Es, fyd, Eyd, DEF, Nd, Mdx, Mdy, TOL_F, TOL_J, TOL_DEF, TOL_k, AREA, Sx, Sy, Ixx, Iyy, Ixy, SINAL_DA_CIRCUICAO, Ast, Tr_I, b, h, Ysmin, Ysmax);
    
    if RUINA_SECAO == false
        DIMENSIONAVEL = true;
        return;
    end

    
    Ast = As_right;
    phi = sqrt(Ast*4/(pi()*Ns));
    As = pi()*phi^2/4*ones(size(Xs));
    [RUINA_SECAO, DEF, NORMA_F] = VerificacaoSecao(Xc, Yc, INC, SIGMAcd, Ec2, Ecu, n, Xs, Ys, As, Nc, Ns, classe_aco, fyk, gamma_s, Es, fyd, Eyd, DEF, Nd, Mdx, Mdy, TOL_F, TOL_J, TOL_DEF, TOL_k, AREA, Sx, Sy, Ixx, Iyy, Ixy, SINAL_DA_CIRCUICAO, Ast, Tr_I, b, h, Ysmin, Ysmax);

    if RUINA_SECAO == true
        fprintf('NAO FOI POSSIVEL DIMENSIONAR A ARMADURA:\n');
        fprintf('AUMENTAR O NUMERO DE BARRAS OU AUMENTAR AS DIMENSOES DA SECAO OU AUMENTAR O fck DO CONCRETO\n');
        return;
    end

    % Inicio da busca pelo metodo da biseccao

    while (I < It_max) && (abs(As_right-As_left)/AREA > TOL_As)
        I = I + 1;
        Ast = psi*As_left + (1-psi)*As_right;
        phi = sqrt(Ast*4/(pi()*Ns));
        As = pi()*phi^2/4*ones(size(Xs));
        [RUINA_SECAO, DEF, NORMA_F] = VerificacaoSecao(Xc, Yc, INC, SIGMAcd, Ec2, Ecu, n, Xs, Ys, As, Nc, Ns, classe_aco, fyk, gamma_s, Es, fyd, Eyd, DEF, Nd, Mdx, Mdy, TOL_F, TOL_J, TOL_DEF, TOL_k, AREA, Sx, Sy, Ixx, Iyy, Ixy, SINAL_DA_CIRCUICAO, Ast, Tr_I, b, h, Ysmin, Ysmax);
        if (RUINA_SECAO == false)
            As_right = Ast;
        else
            As_left = Ast;
        end
    end

    Ast = As_right;
    phi = sqrt(Ast*4/(pi()*Ns));
    As = pi()*phi^2/4*ones(size(Xs));
    for I = 1:5
        fprintf('\n');
    end
    fprintf('DIMENSIONAMENTO DA AREA DE ARMADURA TOTAL PARA O ARRANJO PRÉ-FIXADO:\n');
    fprintf('\t Ast: %.2f mm²\n', 100*sum(As));

end