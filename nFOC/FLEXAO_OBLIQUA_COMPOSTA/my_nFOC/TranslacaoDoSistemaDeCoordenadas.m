function [Xc, Yc, Xs, Ys] = TranslacaoDoSistemaDeCoordenadas(Xc, Yc, Xs, Ys)
    % Entradas:
    % Xc: Coordenadas x dos vértices do contorno
    % Yc: Coordenadas y dos vértices do contorno
    % Xs: Coordenadas x dos pontos do aço
    % Ys: Coordenadas y dos pontos do aço
    
    % Cálculo dos parâmetros da seção
    Nc = size(Xc, 1);
    Ns = size(Xs, 1);
    
    AREA = 0;
    Sx = 0;
    Sy = 0;
    
    for I = 1:Nc
        ai = Xc(I, 1) * Yc(I, 2) - Xc(I, 2) * Yc(I, 1);
        AREA = AREA + ai;
        Sx = Sx + ai * (Yc(I, 1) + Yc(I, 2));
        Sy = Sy + ai * (Xc(I, 1) + Xc(I, 2));
    end

    AREA = (1 / 2) * AREA;
    Sx = (1 / 6) * Sx;
    Sy = (1 / 6) * Sy;
    
    Xcg = Sy / AREA;
    Ycg = Sx / AREA;

    % Translação das coordenadas
    for I = 1:Nc
        Xc(I, 1) = Xc(I, 1) - Xcg;
        Xc(I, 2) = Xc(I, 2) - Xcg;
        Yc(I, 1) = Yc(I, 1) - Ycg;
        Yc(I, 2) = Yc(I, 2) - Ycg;
    end

    for I = 1:Ns
        Xs(I) = Xs(I) - Xcg;
        Ys(I) = Ys(I) - Ycg;
    end
end
