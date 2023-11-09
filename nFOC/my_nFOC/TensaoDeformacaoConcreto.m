function SIGMA = TensaoDeformacaoConcreto(E, SIGMAcd, Ec2, n)

% Initialize SIGMA
SIGMA = 0.0;

% Check if E is less than 0
if E < 0
    SIGMA = 0.0;
end

% Check if E is between 0 and Ec2
if E >= 0 && E < Ec2
    SIGMA = SIGMAcd * (1 - (1 - E / Ec2)^n);
end

% Check if E is greater than or equal to Ec2
if E >= Ec2
    SIGMA = SIGMAcd;
end

end
