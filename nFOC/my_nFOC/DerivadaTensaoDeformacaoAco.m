function DSIGMADE = DerivadaTensaoDeformacaoAco(E, classe_aco, Eyd, fyd)

% Initialize DSIGMADE
DSIGMADE = 0.0;

% Determine behavior based on 'classe_aco'
switch classe_aco
    case 'A'
        if abs(E) <= Eyd
            DSIGMADE = 1 / Eyd;
        end
    case 'B'
        if abs(E) <= 0.7 * Eyd
            DSIGMADE = 1 / Eyd;
        elseif abs(E) > (Eyd + 2)
            DSIGMADE = 0.0;
        else
            DSIGMADE = 3 / sqrt(800 * abs(E) + Eyd * (9 * Eyd - 560));
        end
end

DSIGMADE = fyd * DSIGMADE;

end
