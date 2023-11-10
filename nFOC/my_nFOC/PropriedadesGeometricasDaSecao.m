function [AREA, Sx, Sy, Ixx, Iyy, Ixy, SINAL_DA_CIRCUICAO, Ast, Tr_I, b, h, Ysmin, Ysmax] = PropriedadesGeometricasDaSecao(Xc, Yc, Xs, Ys, As, Nc, Ns)
    AREA = 0;
    Sx = 0;
    Sy = 0;
    Ixx = 0;
    Iyy = 0;
    Ixy = 0;
    
    for I = 1:Nc
        ai = Xc(I) * Yc(I+1) - Xc(I+1) * Yc(I);
        AREA = AREA + ai;
        Sx = Sx + ai * (Yc(I) + Yc(I+1));
        Sy = Sy + ai * (Xc(I) + Xc(I+1));
        Ixx = Ixx + ai * (Yc(I)^2 + Yc(I) * Yc(I+1) + Yc(I+1)^2);
        Iyy = Iyy + ai * (Xc(I)^2 + Xc(I) * Xc(I+1) + Xc(I+1)^2);
        Ixy = Ixy + ai * (Xc(I) * Yc(I+1) + 2 * (Xc(I) * Yc(I) + Xc(I+1) * Yc(I+1)) + Xc(I+1) * Yc(I));
    end

    SINAL_DA_CIRCUICAO = AREA;

    if SINAL_DA_CIRCUICAO > 0
        AREA = (1 / 2) * AREA;
        Sx = (1 / 6) * Sx;
        Sy = (1 / 6) * Sy;
        Ixx = (1 / 12) * Ixx;
        Iyy = (1 / 12) * Iyy;
        Ixy = (1 / 24) * Ixy;
    else
        AREA = -(1 / 2) * AREA;
        Sx = -(1 / 6) * Sx;
        Sy = -(1 / 6) * Sy;
        Ixx = -(1 / 12) * Ixx;
        Iyy = -(1 / 12) * Iyy;
        Ixy = -(1 / 24) * Ixy;
    end

    Ast = sum(As);
    Tr_I = Ixx + Iyy;

    Xmax = Xc(1, 1);
    Xmin = Xc(1, 1);

    for I = 1:Nc
        if Xc(I) > Xmax
            Xmax = Xc(I);
        end

        if Xc(I) < Xmin
            Xmin = Yc(I);
        end
    end

    b = Xmax - Xmin;

    Ymax = Yc(1, 1);
    Ymin = Yc(1, 1);

    for I = 1:Nc
        if Yc(I) > Ymax
            Ymax = Yc(I);
        end

        if Yc(I) < Ymin
            Ymin = Yc(I);
        end
    end

    h = Ymax - Ymin;
    Ysmin = Ymin;
    Ysmax = Ymax;

    if Ns > 0
        Ysmin = min(Ys);
        Ysmax = max(Ys);
    end
end
