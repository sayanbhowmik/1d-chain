function [S] = poissonSolve(S)
% Solve A*phi = RHS

t1 = tic;

RHS = (4*pi/S.epsilon_elec)*(S.b+S.rho);
A = -S.Lap_std + (S.kappa^2 * speye(S.N));

S.phi = A \ RHS;

% S.phi = S.phi - dot(S.W,S.phi)/sum(S.W); % Make it zero mean

fprintf(' Poisson problem took %fs\n',toc(t1));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Using the kernel
% function [S] = poissonSolve(S)
% % Solve A*phi = RHS
% 
% t1 = tic;
% 
% % Calculate Yukawa Kernel
% Kxy = yukawa_kernel(S.x,S.x,S);
% 
% S.phi = Kxy*(S.W.*(S.b+S.rho));
% % S.phi = S.phi - dot(S.W,S.phi)/sum(S.W); % Make it zero mean
% 
% fprintf(' Kernel method took %fs\n',toc(t1));
% end

% function K = yukawa_kernel(x, y, S)
% % Takes two column vectors x and y and returns matrix K(i,j) = K(x_i, y_j)
% 
% % Ensure x and y are column vectors
% x = x(:);
% y = y(:);
% 
% % Calculate distance matrix |x_i - y_j|
% % x is (Nx1), y' is (1xN). Matlab expands this to (NxN)
% raw_dist = abs(x - y');
% dist_matrix = min(raw_dist, S.L - raw_dist);
% 
% % Calculate Kernel
% prefactor = (2 * pi) / (S.kappa * S.epsilon_elec);
% K = prefactor * exp(-S.kappa * dist_matrix);
% end

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