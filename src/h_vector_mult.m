function Hx = h_vector_mult(S,X)
Hx = zeros(size(X));
Veff = S.Veff;

% Kinetic energy laplacian
Laplacian = S.Lap_std;

Hx = Hx - 0.5*Laplacian*X;

% apply Veff
Hx = Hx + Veff.*X;
end