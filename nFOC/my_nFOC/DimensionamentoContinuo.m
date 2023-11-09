function DimensionamentoContinuo()
    %% ?? Precisa passar como argumento as variaveis

    % Define variáveis iniciais
    Nite_DIMENSIONAMENTO = 10000;
    DIMENSIONAVEL = false;

    % Verifica se há barras na seção
    if Ns ~= 0
        Asd = (4 / 100) * AREA;

        if Nds > 0
            if 0.15 * Nds / fyd > (0.4 / 100) * AREA
                Ase = 0.15 * Nds / fyd;
            else
                Ase = (0.4 / 100) * AREA;
            end
        else
            Ase = 0;
        end

        for I = 1:Nite_DIMENSIONAMENTO
            As = Ase / Ns;
            VerificacaoSecao(Nds, Mdxs, Mdys);

            if ~RUINA_SECAO && ~ELU
                As = Ase / Ns;
                DIMENSIONAVEL = true;
                break;
            else
                As = Asd / Ns;
                VerificacaoSecao(Nds, Mdxs, Mdys);

                if ~RUINA_SECAO && ~ELU
                    DIMENSIONAVEL = true;

                    if (Asd - Ase) <= PRECISAO_As
                        As = Asd / Ns;
                        break;
                    else
                        As = (Ase + Asd) / 2 / Ns;
                        VerificacaoSecao(Nds, Mdxs, Mdys);

                        if ~RUINA_SECAO && ~ELU
                            Asd = As * Ns;
                        else
                            Ase = As * Ns;
                        end
                    end
                else
                    DIMENSIONAVEL = false;
                    break;
                end
            end
        end

        Ast = 0;

        for I = 1:NS
            Ast = Ast + As(I);
        end

        Omega = Ast * fyd / (SIGMAcd * AREA);

        if DIMENSIONAVEL
            fprintf('DIMENSIONAMENTO DA AREA DE ARMADURA TOTAL PARA O ARRANJO PRÉ-FIXADO:\n');
            fprintf('   As: %.2f cm²\n', Ast);
        else
            fprintf('NAO FOI POSSIVEL DIMENSIONAR A ARMADURA:\n');
            fprintf('AUMENTAR O NUMERO DE BARRAS OU AUMENTAR AS DIMENSOES DA SECAO OU AUMENTAR O fck DO CONCRETO\n');
        end
    else
        fprintf('O DIMENSIONAMENTO REQUER PELO MENOS UMA BARRA PARA O ARRANJO PRÉ-FIXADO:\n');
    end
end
