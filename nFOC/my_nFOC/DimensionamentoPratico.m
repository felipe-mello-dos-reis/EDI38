function [DIMENSIONAVEL, phi, Ast, DEF, NORMA_F] = DimensionamentoPratico(Xc, Yc, INC, SIGMAcd, Ec2, Ecu, n, Xs, Ys, As, Nc, Ns, classe_aco, fyk, gamma_s, Es, fyd, Eyd, DEF, Nd, Mdx, Mdy, TOL_F, TOL_J, TOL_DEF, TOL_k, AREA, Sx, Sy, Ixx, Iyy, Ixy, SINAL_DA_CIRCUICAO, Ast, Tr_I, b, h, Ysmin, Ysmax)

    DIMENSIONAVEL = false;
    phis = [0.5, 0.63, 0.8, 1.0, 1.25, 1.6, 2.0, 2.5, 3.2, 4.0];
    for I = 1:length(phis)
        phi = phis(I);
        fprintf("-> phi = %.1f mm\n", phi*10);
        As = pi()*phi^2/4*ones(size(Xs));
        Ast = sum(As);
        [RUINA_SECAO, DEF, NORMA_F] = VerificacaoSecao(Xc, Yc, INC, SIGMAcd, Ec2, Ecu, n, Xs, Ys, As, Nc, Ns, classe_aco, fyk, gamma_s, Es, fyd, Eyd, DEF, Nd, Mdx, Mdy, TOL_F, TOL_J, TOL_DEF, TOL_k, AREA, Sx, Sy, Ixx, Iyy, Ixy, SINAL_DA_CIRCUICAO, Ast, Tr_I, b, h, Ysmin, Ysmax);
        if RUINA_SECAO == false
            DIMENSIONAVEL = true;
            break
        end
    end
    
    for I = 1:10
        fprintf('\n');
    end

    if DIMENSIONAVEL == true
        fprintf('DIMENSIONAMENTO DA AREA DE ARMADURA TOTAL PARA O ARRANJO PRÉ-FIXADO:\n');
        fprintf('\t Ast: %.2f mm²\n', 100*sum(As));
        fprintf("\t phi = %.1f mm\n", phi*10);
    else
        fprintf('NAO FOI POSSIVEL DIMENSIONAR A ARMADURA:\n');
        fprintf('AUMENTAR O NUMERO DE BARRAS OU AUMENTAR AS DIMENSOES DA SECAO OU AUMENTAR O fck DO CONCRETO\n');
    end
end