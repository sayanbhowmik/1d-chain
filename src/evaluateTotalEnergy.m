function [Etot,Eband,Exc,Exc_dc,Eelec_dc,Eent] = evaluateTotalEnergy(S)

% Band structure energy
Eband = 0;
for spin = 1:S.nspin
    Eband = Eband + S.occfac * sum(S.EigVal(:,1).*S.occ(:,1)) ;
end

% Exchange-correlation energy
rho = S.rho;
% Check if density is too small
INDX_zerorho = (rho < S.xc_rhotol);
rho(INDX_zerorho) = S.xc_rhotol;

Exc = sum(S.e_xc.*(rho).*S.W);
Exc = S.XCswitch*Exc;
% Exchange-correlation energy double counting correction
Exc_dc = sum(S.Vxc.*rho.*S.W) ;
Exc_dc = S.XCswitch*Exc_dc;

% Electrostatic energy double counting correction
Eelec_dc = 0.5*sum((S.b-S.rho(:,1)).*S.phi.*S.W);

% Electronic entropy
Eent = 0 ;
ks = 1;
for spin = 1:S.nspin
    % fermi-dirac smearing
    Eent_v = S.occfac*(1/S.bet)*(S.occ(:,ks).*log(S.occ(:,ks))+(1-S.occ(:,ks)).*log(1-S.occ(:,ks)));
    Eent_v(isnan(Eent_v)) = 0.0 ;
    Eent = Eent + sum(Eent_v);
end

% Total free energy
Etot = Eband + Exc - Exc_dc + Eelec_dc + Eent - S.Eself;

fprintf(2,' ------------------\n');
fprintf(' Eband = %.8f\n', Eband);
fprintf(' Exc = %.8f\n', Exc);
fprintf(' Exc_dc = %.8f\n', Exc_dc);
fprintf(' Eelec_dc = %.8f\n', Eelec_dc);
fprintf(' Eent = %.8f\n', Eent);
% fprintf(' E_corr = %.8f\n', S.E_corr);
fprintf(' Eself = %.8f\n', S.Eself);
fprintf(' Etot = %.8f\n', Etot);
fprintf(2,' ------------------\n');
end