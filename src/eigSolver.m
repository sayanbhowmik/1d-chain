function [psi, EigVal] = eigSolver(S)
% eigen problem size
N = S.N*S.nspinor_eig;

rng('default'); % Initialize random number generator
opts = struct('maxit', 1000, 'tol', S.TOL_LANCZOS, 'v0', rand(N,1),'issym', 1);

Hfun = @(x) h_vector_mult(S,x);

[psi, EigVal] = eigs(Hfun,N,S.Nev,'sr',opts);
EigVal = diag(EigVal);

% Sort eigenvalues in ascending order
[EigVal, idx] = sort(EigVal, 'ascend');

% Reorder the eigenvectors to match the sorted eigenvalues
psi = psi(:, idx);

% Normalize psi  s.t. integral(psi' * psi) = 1
psi = psi / sqrt(S.dx);

end