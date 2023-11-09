function J = DerivadasDosEsforcosResistentes(DEF, Xc, Yc, INC, Nc, Xs, Ys, As, Ns, TOL_k, TOL_DEF, SINAL_DA_CIRCUICAO, AREA, Sx, Sy, Ixx, Iyy, Ixy, ERc)

%% ?? Verifique f10i

% Initialize J matrices
Jc = zeros(3, 3);
Js = zeros(3, 3);

INPUTS;

% Loop for concrete contribution
for I = 1:Nc
    Dxi = Xc(INC(I, 2)) - Xc(INC(I, 1));
    Dyi = Yc(INC(I, 2)) - Yc(INC(I, 1));
    Ei = DEF(1) + DEF(3) * Xc(INC(I, 1)) - DEF(2) * Yc(INC(I, 1));
    Ei_1 = DEF(1) + DEF(3) * Xc(INC(I, 2)) - DEF(2) * Yc(INC(I, 2));
    DEi = Ei_1 - Ei;
    hi = Xc(INC(I, 1)) * Ei_1 - Xc(INC(I, 2)) * Ei;
    gi = Yc(INC(I, 1)) * Ei_1 - Yc(INC(I, 2)) * Ei;
    
    [I1i, I2i, I3i, J1i, J2i, K1i] = IntegraisDiagramaTensaoDeformacao(Ei);
    [I1i_1, I2i_1, I3i_1, J1i_1, J2i_1, K1i_1] = IntegraisDiagramaTensaoDeformacao(Ei_1);
    
    [SIGMAi, DSIGMADEi] = TensaoDeformacaoConcreto(Ei, classe_aco, Eyd, fyd);

    DI1i = I1i_1 - I1i;
    DJ1i = J1i_1 - J1i;
    DJ2i = J2i_1 - J2i;
    
    if abs(DEi) <= TOL_DEF
        f6i = SIGMAi;
        f7i = SIGMAi * (Yc(INC(I, 1)) + Yc(INC(I, 2))) / 2;
        f8i = SIGMAi * (Xc(INC(I, 1)) + Xc(INC(I, 2))) / 2;
        f9i = SIGMAi * (Yc(INC(I, 1))^2 + Yc(INC(I, 1)) * Yc(INC(I, 2)) + Yc(INC(I, 2))^2) / 3;
        % ?? ACHO QUE TEM PARENTESES A MAIS AQUI
        f10i = SIGMAi * (Xc(INC(I, 1)) * Yc(INC(I, 2)) + 2 * (Xc(INC(I, 1)) * Yc(INC(I, 1)) + Xc(INC(I, 2)) * Yc(INC(I, 2)) + Xc(INC(I, 2)) * Yc(INC(I, 1)))) / 6;
        f11i = SIGMAi * (Xc(INC(I, 1))^2 + Xc(INC(I, 1)) * Xc(INC(I, 2)) + Xc(INC(I, 2))^2) / 3;
    else
        f6i = DI1i / DEi;
        f7i = (gi * DI1i + Dyi * DJ1i) / DEi^2;
        f8i = (hi * DI1i + Dxi * DJ1i) / DEi^2;
        f9i = (gi^2 * DI1i + 2 * gi * Dyi * DJ1i + Dyi^2 * DJ2i) / DEi^3;
        f10i = (hi * gi * DI1i + (gi * Dxi + hi * Dyi) * DJ1i + Dxi * Dyi * DJ2i) / DEi^3;
        f11i = (hi^2 * DI1i + 2 * hi * Dxi * DJ1i + Dxi^2 * DJ2i) / DEi^3;
    end
    
    if SINAL_DA_CIRCUICAO > 0
        Jc(1, 1) = Jc(1, 1) + Dxi * f6i;
        Jc(1, 2) = Jc(1, 2) - Dxi * f7i;
        Jc(1, 3) = Jc(1, 3) + Dxi * f8i;
        Jc(2, 1) = Jc(2, 1) - Dxi * f7i;
        Jc(2, 2) = Jc(2, 2) + Dxi * f9i;
        Jc(2, 3) = Jc(2, 3) - Dxi * f10i;
        Jc(3, 1) = Jc(3, 1) + Dxi * f8i;
        Jc(3, 2) = Jc(3, 2) - Dxi * f10i;
        Jc(3, 3) = Jc(3, 3) + Dxi * f11i;
    else
        Jc(1, 1) = Jc(1, 1) - Dxi * f6i;
        Jc(1, 2) = Jc(1, 2) + Dxi * f7i;
        Jc(1, 3) = Jc(1, 3) - Dxi * f8i;
        Jc(2, 1) = Jc(2, 1) + Dxi * f7i;
        Jc(2, 2) = Jc(2, 2) - Dxi * f9i;
        Jc(2, 3) = Jc(2, 3) + Dxi * f10i;
        Jc(3, 1) = Jc(3, 1) - Dxi * f8i;
        Jc(3, 2) = Jc(3, 2) + Dxi * f10i;
        Jc(3, 3) = Jc(3, 3) - Dxi * f11i;
    end
end

% Update Jc matrix
Jc = Jc * (1 / DEF(2));

% Loop for steel contribution
for I = 1:Ns
    Ei = DEF(1) + DEF(3) * Xs(I) - DEF(2) * Ys(I);
    DSIGMADEi = DerivadaTensaoDeformacaoAco(Ei, classe_aco, Eyd, fyd);

    Js(1, 1) = Js(1, 1) + DSIGMADEi * As(I);
    Js(1, 2) = Js(1, 2) - DSIGMADEi * Ys(I) * As(I);
    Js(1, 3) = Js(1, 3) + DSIGMADEi * Xs(I) * As(I);
    Js(2, 1) = Js(2, 1) - DSIGMADEi * Ys(I) * As(I);
    Js(2, 2) = Js(2, 2) + DSIGMADEi * (Ys(I)^2) * As(I);
    Js(2, 3) = Js(2, 3) - DSIGMADEi * Xs(I) * Ys(I) * As(I);
    Js(3, 1) = Js(3, 1) + DSIGMADEi * Xs(I) * As(I);
    Js(3, 2) = Js(3, 2) - DSIGMADEi * Xs(I) * Ys(I) * As(I);
    Js(3, 3) = Js(3, 3) + DSIGMADEi * (Xs(I)^2) * As(I);
end

% Update Js matrix
J = Jc + Js;

end
