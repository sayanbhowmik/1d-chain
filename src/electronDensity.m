function S = electronDensity(S)

rho_d = zeros(S.N,S.nspinor);
for spinor = 1:S.nspinor
    ndrange = (1+(spinor-1)*S.N:spinor*S.N);
    rho_d(:,spinor) = rho_d(:,spinor) + S.occfac* sum( S.psi(ndrange,:).*conj(S.psi(ndrange,:)).*S.occ(:,1)',2);
end
rho_d = real(rho_d);
rhotot = sum(rho_d,2);
S.rho = rhotot;
end