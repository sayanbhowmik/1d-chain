function S = calculate_Eself(S)
S.Eself = 0;

% Calculate Yukawa Kernel
Kxy = yukawa_kernel(S.x,S.x,S);

atom_positions = S.Atoms;
for JJ_a = 1 : S.n_atm
    % Calculate Signed Minimum Image Distance
    dx_raw = S.x - atom_positions(JJ_a);
    % This trick automatically wraps the distance to [-L/2, L/2]
    dx_pbc = dx_raw - S.L * round(dx_raw / S.L);

    % Calculate Gaussian
    gauss = -(S.Z(JJ_a) / sqrt(2*pi*S.b_sigma(JJ_a)^2)) * exp(-dx_pbc.^2 / (2*S.b_sigma(JJ_a)^2));

    % Evaluate: \int \int K(x,y) b_{JJ_a}(x) b_{JJ_a}(y) dx dy
    phi_self = Kxy*(S.W.*gauss);
    integral = sum(S.W.*phi_self.*gauss);

    S.Eself = S.Eself + 0.5*integral;
end

% S.Eself = 0.5*S.Eself;

end

function K = yukawa_kernel(x, y, S)
% Takes two column vectors x and y and returns matrix K(i,j) = K(x_i, y_j)

% Ensure x and y are column vectors
x = x(:);
y = y(:);

% Calculate distance matrix |x_i - y_j|
% x is (Nx1), y' is (1xN). Matlab expands this to (NxN)
raw_dist = abs(x - y');
dist_matrix = min(raw_dist, S.L - raw_dist);

% Calculate Kernel
prefactor = (2 * pi) / (S.kappa * S.epsilon_elec);
K = prefactor * exp(-S.kappa * dist_matrix);
end

% function K = yukawa_kernel(x,y,S)
% % Infinite periodic case using cosh and sinh
% 
% % Takes two vectors x and y and returns kernel K(x,y)
% % Distance matrix with Minimum Image Convention
% x = x(:);
% y = y(:);
% raw_dist = abs(x - y');
% dist_matrix = min(raw_dist, S.L - raw_dist);
% 
% % Periodic Kernel Formula (Cosh)
% % This sums infinite periodic images exactly
% numerator = cosh(S.kappa * (0.5 * S.L - dist_matrix));
% denominator = 2 * S.kappa * sinh(S.kappa * S.L / 2);
% 
% prefactor = 4 * pi / S.epsilon_elec; % Note prefactor change for this form
% K = prefactor * (numerator / denominator);
% end