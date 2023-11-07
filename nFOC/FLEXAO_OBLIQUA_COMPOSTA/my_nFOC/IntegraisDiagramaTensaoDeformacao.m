function [I1, I2, I3, J1, J2, K1] = IntegraisDiagramaTensaoDeformacao(E, SIGMAcd, n, Ec2)
    %% ?? Verificar I3 e K1

    t = zeros(1, 3);

    if (E >= 0 && E < Ec2)
        t(1) = (1 / (n + 1)) * (1 - (1 - E / Ec2)^(n + 1));
        t(2) = (1 / (n + 2)) * (1 - (1 - E / Ec2)^(n + 2));
        t(3) = (1 / (n + 3)) * (1 - (1 - E / Ec2)^(n + 3));
    end

    if (E >= Ec2)
        t(1) = 1 / (n + 1);
        t(2) = 1 / (n + 2);
        t(3) = 1 / (n + 3);
    end

    if (E < 0)
        I1 = SIGMAcd * 0;
    end

    if (E >= 0)
        I1 = SIGMAcd * (E - Ec2 * t(1));
    end

    if (E < 0)
        I2 = SIGMAcd * 0;
    end

    if (E >= 0)
        I2 = SIGMAcd * (E^2 / 2 - (Ec2 / (n + 1)) * (E - Ec2 * t(2)));
    end

    if (E < 0)
        I3 = SIGMAcd * 0;
    end

    if (E >= 0)
        %% ?? parenteses a mais
        I3 = SIGMAcd * (E^3 / 6 - (Ec2 / (n + 1)) * (E^2 / 2 - (Ec2 / (n + 2)) * (E - Ec2 * t(3))));
    end

    if (E < 0)
        J1 = SIGMAcd * 0;
    end

    if (E >= 0)
        J1 = SIGMAcd * (E^2 / 2 - Ec2^2 * (t(1) - t(2)));
    end

    if (E < 0)
        J2 = SIGMAcd * 0;
    end

    if (E >= 0)
        J2 = SIGMAcd * (E^3 / 3 - Ec2^3 * (t(1) - 2 * t(2) + t(3)));
    end

    if (E < 0)
        K1 = SIGMAcd * 0;
    end

    if (E >= 0)
        %% ?? parenteses a mais
        K1 = SIGMAcd * (E^3 / 3 - (Ec2 / (n + 1)) * (E^2 / 2 - Ec2^2 * (t(2) - t(3))));
    end
end
