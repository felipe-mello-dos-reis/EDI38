function J = DerivadasDosEsforcosResistentes(Xc, Yc, INC, SIGMAcd, Ec2, Ecu, n, Xs, Ys, As, Nc, Ns, classe_aco, fyk, gamma_s, Es, fyd, Eyd, DEF, Nd, Mdx, Mdy, TOL_F, TOL_J, TOL_DEF, TOL_k, AREA, Sx, Sy, Ixx, Iyy, Ixy, SINAL_DA_CIRCUICAO, Ast, Tr_I, b, h, Ysmin, Ysmax)
% Calculando os esforços resistentes
[ER, ERc, ERs] = EsforcosResistentes(Xc, Yc, INC, SIGMAcd, Ec2, Ecu, n, Xs, Ys, As, Nc, Ns, classe_aco, fyk, gamma_s, Es, fyd, Eyd, DEF, Nd, Mdx, Mdy, TOL_F, TOL_J, TOL_DEF, TOL_k, AREA, Sx, Sy, Ixx, Iyy, Ixy, SINAL_DA_CIRCUICAO, Ast, Tr_I, b, h, Ysmin, Ysmax);


% Derivadas dos esforços resistentes: contribuição do concreto
if abs(DEF(2)*h) >= TOL_k
    Jc = zeros(3);
    % Js = zeros(3);
    Jc(1, 2) = -ERc(1);
    Jc(2, 1) = -ERc(1);
    Jc(2, 2) = -2 * ERc(2);
    Jc(2, 3) = -ERc(3);
    Jc(3, 2) = -ERc(3);
    
    for I = 1:Nc
        Dxi = Xc(INC(I, 2)) - Xc(INC(I, 1));
        Dyi = Yc(INC(I, 2)) - Yc(INC(I, 1));
        Ei = DEF(1) + DEF(3)*Xc(INC(I, 1)) - DEF(2)*Yc(INC(I, 1));
        Ei_1 = DEF(1) + DEF(3)*Xc(INC(I, 2)) - DEF(2)*Yc(INC(I, 2));
        DEi = Ei_1 - Ei;
        hi = Xc(INC(I, 1))*Ei_1 - Xc(INC(I, 2))*Ei;
        gi = Yc(INC(I, 1))*Ei_1 - Yc(INC(I, 2))*Ei;
        
        % IntegraisDiagramaTensaoDeformacao(Ei, I1i, I2i, I3i, J1i, J2i, K1i);
        [I1i, I2i, I3i, J1i, J2i, K1i] = IntegraisDiagramaTensaoDeformacao(Ei, SIGMAcd, n, Ec2);
        % IntegraisDiagramaTensaoDeformacao(Ei_1, I1i_1, I2i_1, I3i_1, J1i_1, J2i_1, K1i_1);
        [I1i_1, I2i_1, I3i_1, J1i_1, J2i_1, K1i_1] = IntegraisDiagramaTensaoDeformacao(Ei_1, SIGMAcd, n, Ec2);
        % TensaoDeformacaoConcreto(Ei, SIGMAi);
        SIGMAi = TensaoDeformacaoConcreto(Ei, SIGMAcd, Ec2, n);
        
        DI1i = I1i_1 - I1i;
        DJ1i = J1i_1 - J1i;
        DJ2i = J2i_1 - J2i;
        
        if abs(DEi) <= TOL_DEF
            f6i = SIGMAi;
            f7i = SIGMAi * (Yc(INC(I, 1)) + Yc(INC(I, 2))) / 2;
            f8i = SIGMAi * (Xc(INC(I, 1)) + Xc(INC(I, 2))) / 2;
            f9i = SIGMAi * (Yc(INC(I, 1))^2 + Yc(INC(I, 1))*Yc(INC(I, 2)) + Yc(INC(I, 2))^2) / 3;
            f10i = SIGMAi * (Xc(INC(I, 1))*Yc(INC(I, 2)) + 2 * (Xc(INC(I, 1))*Yc(INC(I, 1)) + Xc(INC(I, 2))*Yc(INC(I, 2))) + Xc(INC(I, 2))*Yc(INC(I, 1))) / 6;
            f11i = SIGMAi * (Xc(INC(I, 1))^2 + Xc(INC(I, 1))*Xc(INC(I, 2)) + Xc(INC(I, 2))^2) / 3;
        else
            f6i = DI1i / DEi;
            f7i = (gi*DI1i + Dyi*DJ1i) / DEi^2;
            f8i = (hi*DI1i + Dxi*DJ1i) / DEi^2;
            f9i = ((gi^2)*DI1i + 2 * gi*Dyi*DJ1i + (Dyi^2)*DJ2i) / DEi^3;
            f10i = (hi*gi*DI1i + (gi*Dxi + hi*Dyi)*DJ1i + Dxi*Dyi*DJ2i) / DEi^3;
            f11i = ((hi^2)*DI1i + 2 * hi*Dxi*DJ1i + (Dxi^2)*DJ2i) / DEi^3;
        end
        
        if SINAL_DA_CIRCUICAO > 0
            Jc(1, 1) = Jc(1, 1) + Dxi*f6i;
            Jc(1, 2) = Jc(1, 2) - Dxi*f7i;
            Jc(1, 3) = Jc(1, 3) + Dxi*f8i;
            Jc(2, 1) = Jc(2, 1) - Dxi*f7i;
            Jc(2, 2) = Jc(2, 2) + Dxi*f9i;
            Jc(2, 3) = Jc(2, 3) - Dxi*f10i;
            Jc(3, 1) = Jc(3, 1) + Dxi*f8i;
            Jc(3, 2) = Jc(3, 2) - Dxi*f10i;
            Jc(3, 3) = Jc(3, 3) + Dxi*f11i;
        else
            Jc(1, 1) = Jc(1, 1) - Dxi*f6i;
            Jc(1, 2) = Jc(1, 2) + Dxi*f7i;
            Jc(1, 3) = Jc(1, 3) - Dxi*f8i;
            Jc(2, 1) = Jc(2, 1) + Dxi*f7i;
            Jc(2, 2) = Jc(2, 2) - Dxi*f9i;
            Jc(2, 3) = Jc(2, 3) + Dxi*f10i;
            Jc(3, 1) = Jc(3, 1) - Dxi*f8i;
            Jc(3, 2) = Jc(3, 2) + Dxi*f10i;
            Jc(3, 3) = Jc(3, 3) - Dxi*f11i;
        end
    end
    Jc = Jc * (1/DEF(2));
elseif abs(DEF(3)*h) >= TOL_k
    Jc = zeros(3);
    Jc(1, 3) = -ERc(1);
    Jc(3, 1) = -ERc(1);
    Jc(2, 3) = -ERc(2);
    Jc(3, 2) = -ERc(2);
    Jc(3, 3) = -2 * ERc(3);
    
    for I = 1:Nc
        Dxi = Xc(INC(I, 2)) - Xc(INC(I, 1));
        Dyi = Yc(INC(I, 2)) - Yc(INC(I, 1));
        Ei = DEF(1) + DEF(3)*Xc(INC(I, 1)) - DEF(2)*Yc(INC(I, 1));
        Ei_1 = DEF(1) + DEF(3)*Xc(INC(I, 2)) - DEF(2)*Yc(INC(I, 2));
        DEi = Ei_1 - Ei;
        hi = Xc(INC(I, 1))*Ei_1 - Xc(INC(I, 2))*Ei;
        gi = Yc(INC(I, 1))*Ei_1 - Yc(INC(I, 2))*Ei;
        
        IntegraisDiagramaTensaoDeformacao(Ei, I1i, I2i, I3i, J1i, J2i, K1i);
        IntegraisDiagramaTensaoDeformacao(Ei_1, I1i_1, I2i_1, I3i_1, J1i_1, J2i_1, K1i_1);
        TensaoDeformacaoConcreto(Ei, SIGMAi);
        
        DI1i = I1i_1 - I1i;
        DJ1i = J1i_1 - J1i;
        DJ2i = J2i_1 - J2i;
        
        if abs(DEi) <= TOL_DEF
            f6i = SIGMAi;
            f7i = SIGMAi * (Yc(INC(I, 1)) + Yc(INC(I, 2))) / 2;
            f8i = SIGMAi * (Xc(INC(I, 1)) + Xc(INC(I, 2))) / 2;
            f9i = SIGMAi * (Yc(INC(I, 1))^2 + Yc(INC(I, 1))*Yc(INC(I, 2)) + Yc(INC(I, 2))^2) / 3;
            f10i = SIGMAi * (Xc(INC(I, 1))*Yc(INC(I, 2)) + 2 * (Xc(INC(I, 1))*Yc(INC(I, 1)) + Xc(INC(I, 2))*Yc(INC(I, 2))) + Xc(INC(I, 2))*Yc(INC(I, 1))) / 6;
            f11i = SIGMAi * (Xc(INC(I, 1))^2 + Xc(INC(I, 1))*Xc(INC(I, 2)) + Xc(INC(I, 2))^2) / 3;
        else
            f6i = DI1i / DEi;
            f7i = (gi*DI1i + Dyi*DJ1i) / DEi^2;
            f8i = (hi*DI1i + Dxi*DJ1i) / DEi^2;
            f9i = ((gi^2)*DI1i + 2 * gi*Dyi*DJ1i + (Dyi^2)*DJ2i) / DEi^3;
            f10i = (hi*gi*DI1i + (gi*Dxi + hi*Dyi)*DJ1i + Dxi*Dyi*DJ2i) / DEi^3;
            f11i = ((hi^2)*DI1i + 2 * hi*Dxi*DJ1i + (Dxi^2)*DJ2i) / DEi^3;
        end
        
        if SINAL_DA_CIRCUICAO > 0
            Jc(1, 1) = Jc(1, 1) + Dyi*f6i;
            Jc(1, 2) = Jc(1, 2) - Dyi*f7i;
            Jc(1, 3) = Jc(1, 3) + Dyi*f8i;
            Jc(2, 1) = Jc(2, 1) - Dyi*f7i;
            Jc(2, 2) = Jc(2, 2) + Dyi*f9i;
            Jc(2, 3) = Jc(2, 3) - Dyi*f10i;
            Jc(3, 1) = Jc(3, 1) + Dyi*f8i;
            Jc(3, 2) = Jc(3, 2) - Dyi*f10i;
            Jc(3, 3) = Jc(3, 3) + Dyi*f11i;
        else
            Jc(1, 1) = Jc(1, 1) - Dyi*f6i;
            Jc(1, 2) = Jc(1, 2) + Dyi*f7i;
            Jc(1, 3) = Jc(1, 3) - Dyi*f8i;
            Jc(2, 1) = Jc(2, 1) + Dyi*f7i;
            Jc(2, 2) = Jc(2, 2) - Dyi*f9i;
            Jc(2, 3) = Jc(2, 3) + Dyi*f10i;
            Jc(3, 1) = Jc(3, 1) - Dyi*f8i;
            Jc(3, 2) = Jc(3, 2) + Dyi*f10i;
            Jc(3, 3) = Jc(3, 3) - Dyi*f11i;
        end
    end
    Jc = Jc * (1/DEF(3));
    
    
elseif abs(DEF(2)*h) <= TOL_k && abs(DEF(3)*h) <= TOL_k
    DSIGMADE0 = DerivadaTensaoDeformacaoConcreto(DEF(1), SIGMAcd, Ec2, n);
    
    Jc = zeros(3, 3);
    Jc(1, 1) = DSIGMADE0 * AREA;
    Jc(1, 2) = DSIGMADE0 * Sx;
    Jc(1, 3) = DSIGMADE0 * Sy;
    Jc(2, 1) = DSIGMADE0 * Sx;
    Jc(2, 2) = DSIGMADE0 * Ixx;
    Jc(2, 3) = -DSIGMADE0 * Ixy;
    Jc(3, 1) = DSIGMADE0 * Sy;
    Jc(3, 2) = -DSIGMADE0 * Ixy;
    Jc(3, 3) = DSIGMADE0 * Iyy;
    
end

Js = zeros(3,3);

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

J = Jc + Js;
end


