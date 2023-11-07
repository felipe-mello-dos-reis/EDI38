function DiagramasEsforcosResistentes()
    %% ?? Precisa passar como argumento as variaveis
    % Define os pontos PA, PB, PC, PD, PE e PF
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
    DEF = [0; 0];

    % Inicializa o índice para preencher os resultados
    II = 0;

    % Loop para preencher os diagramas de Nr e Mrx para as arestas AB, BC, CD, DE, EF e FA
    for I = 1:NPONTOS_DIAGRAMA
        % Calcula o ponto intermediário P para a aresta AB
        P = PA + (PB - PA) * I / NPONTOS_DIAGRAMA;
        DEF(1) = P(2);
        DEF(2) = P(1);
        [ER, ~] = EsforcosResistentes(DEF); % Substitua com a função correta
        II = II + 1;
        DIAGRAMA_Nr_Mrx(II, 1) = ER(1);
        DIAGRAMA_Nr_Mrx(II, 2) = ER(2);
    end

    % Repita o loop acima para as outras arestas (BC, CD, DE, EF, FA)

    % Crie um arquivo de saída e escreva os resultados
    fid = fopen('DIAGRAMAS_ESFORCOS_RESISTENTES.txt', 'w');
    fprintf(fid, '******************** DIAGRAMAS DE ESFORCOS RESISTENTES ********************\n');
    fprintf(fid, 'PONTOS DO DIAGRAMA Nr Mrx\n\n');
    fprintf(fid, '   PONTO          Nr               Mrx\n');
    
    for I = 1:6 * NPONTOS_DIAGRAMA
        fprintf(fid, '%5d     %13.6E     %13.6E\n', I, DIAGRAMA_Nr_Mrx(I, 1), DIAGRAMA_Nr_Mrx(I, 2));
    end
    
    fclose(fid);
end
