function DSIGMADE = DerivadaTensaoDeformacaoConcreto(E, SIGMAcd, Ec2, n)

% Check if E is less than 0
if E < 0
    DSIGMADE = 0.0;
end

% Check if E is between 0 and Ec2
if E >= 0 && E < Ec2
    DSIGMADE = SIGMAcd * (n / Ec2) * (1 - E / Ec2)^(n - 1);
end

% Check if E is greater than or equal to Ec2
if E >= Ec2
    DSIGMADE = 0.0;
end

end
