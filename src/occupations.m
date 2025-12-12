function S = occupations(S)

FermiEnergyEvaluator = @(lambda_f_g) fermiCalc(lambda_f_g,S);
S.lambda_f = fzero(FermiEnergyEvaluator,0);
fprintf(' Fermi energy = %f\n',S.lambda_f);

% fermi-dirac smearing
S.occ = 1./(1+exp(S.bet*(S.EigVal-S.lambda_f)));
end

function f = fermiCalc(lambda_f_g, S)
f = 0;
ks = 1;
for spin = 1:S.nspin
    % fermi-dirac smearing
    f = f + S.occfac * sum(1./(1+exp(S.bet*(S.EigVal(:,ks)-lambda_f_g)))) ;
end
f = f + S.NegCharge;
end