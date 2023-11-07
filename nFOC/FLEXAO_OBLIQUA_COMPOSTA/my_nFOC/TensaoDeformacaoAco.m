function SIGMA = TensaoDeformacaoAco(E, classe_aco, Fyd, Eyd)


    % Verifique a classe do aço e atribua os valores apropriados
    switch classe_aco
        case 'A'
            if abs(E) <= Eyd
                % Tensão é calculada como a tensão de escoamento multiplicada pela deformação
                SIGMA = Fyd * (E / Eyd);
            else
                % Lidar com diferentes casos quando E é maior que Eyd (escoamento)
                if E < 0
                    % Tensão é -Fyd
                    SIGMA = -Fyd;
                elseif E == 0
                    % Tensão é 0
                    SIGMA = 0;
                else
                    % Tensão é Fyd
                    SIGMA = Fyd;
                end
            end
        case 'B'
            if abs(E) <= 0.7 * Eyd
                % Tensão é calculada como a tensão de escoamento multiplicada pela deformação
                SIGMA = Fyd * (E / Eyd);
            elseif abs(E) > (Eyd + 2)
                % Lidar com diferentes casos quando E é maior que Eyd (escoamento)
                if E < 0
                    % Tensão é -Fyd
                    SIGMA = -Fyd;
                elseif E == 0
                    % Tensão é 0
                    SIGMA = 0;
                else
                    % Tensão é Fyd
                    SIGMA = Fyd;
                end
            else
                % Lidar com um caso mais complexo para E entre 0.7 * Eyd e (Eyd + 2)
                if E < 0
                    % Use uma fórmula complexa para calcular a tensão
                    SIGMA = -Fyd * (280 - 9 * Eyd + 3 * sqrt(800 * abs(E) + Eyd * (9 * Eyd - 560))) / 400;
                elseif E == 0
                    % Tensão é 0
                    SIGMA = 0;
                else
                    % Use uma fórmula complexa para calcular a tensão
                    SIGMA = Fyd * (280 - 9 * Eyd + 3 * sqrt(800 * abs(E) + Eyd * (9 * Eyd - 560))) / 400;
                end
            end
    end
end
