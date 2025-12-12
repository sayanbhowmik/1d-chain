function S = calculate_Eself(S)
S.Eself = 0;

% Poisson way
for JJ_a = 1 : S.n_atm
    % Initialize b_j on full grid size
    bJ_full = zeros(S.N,1);

    % Atom position of atom JJ_a
    x0 = S.Atoms(JJ_a,1);

    % Find index of rb_x for the atom type
    idx_type = find(S.unique_Z == S.Z(JJ_a));

    % Periodic boundary
    n_image_xl = floor((S.Atoms(JJ_a,1) + S.rb_x(idx_type))/S.L);
    n_image_xr = floor((S.L - S.Atoms(JJ_a,1)+S.rb_x(idx_type)-S.dx)/S.L);

    % Total number of images of atom JJ_a (including atom JJ_a)
    n_image_total = (n_image_xl+n_image_xr+1) ;

    % Find the coordinates for all the images
	xx_img = (-n_image_xl : n_image_xr) * S.L + x0;

    % Loop over all image(s) of atom JJ_a (including atom JJ_a)
    for count_image = 1:n_image_total

        x0_i = xx_img(count_image);

        % Starting and ending indices of b-region
        ii_s = ceil ((x0_i - S.rb_x(idx_type))/S.dx) + 1;
        ii_e = floor((x0_i + S.rb_x(idx_type))/S.dx) + 1;

        ii_s = max(ii_s,1);
        ii_e = min(ii_e,S.N);

        xx = (ii_s-S.FDn-1:ii_e+S.FDn-1)*S.dx;
        dd = bsxfun(@minus,xx,x0_i)';
        % dd = sqrt(dd.^2);

        % Pseudocharge density calculation
		II = 1+S.FDn : size(dd,1)-S.FDn;

        bJ = pseudochargeDensity_atom(dd,II,S.Z(JJ_a),S.b_sigma(JJ_a));

        II_act = ii_s:ii_e;

        bJ_full(II_act) = bJ_full(II_act) + bJ(II);
    end
    RHS = (4*pi/S.epsilon_elec)*(bJ_full);
    A = -S.Lap_std + (S.kappa^2 * speye(S.N));
    
    phi_bj = A \ RHS;

    eself_j = 0.5*sum(S.W.*bJ_full.*phi_bj);
    S.Eself = S.Eself+eself_j;
end

fprintf("Eself = %0.14f\n",S.Eself);

% %%%%%%%%%%%%%%%%%%%%%%%% Yukawa Kernel way %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Likely generates integration errors for cusp condition at |x-y| at x=y
% % Calculate Yukawa Kernel
% Kxy = yukawa_kernel(S.x,S.x,S);
% 
% atom_positions = S.Atoms;
% for JJ_a = 1 : S.n_atm
%     % Calculate Signed Minimum Image Distance
%     dx_raw = S.x - atom_positions(JJ_a);
%     % This trick automatically wraps the distance to [-L/2, L/2]
%     dx_pbc = dx_raw - S.L * round(dx_raw / S.L);
% 
%     % Calculate Gaussian
%     gauss = -(S.Z(JJ_a) / sqrt(2*pi*S.b_sigma(JJ_a)^2)) * exp(-dx_pbc.^2 / (2*S.b_sigma(JJ_a)^2));
% 
%     % Evaluate: \int \int K(x,y) b_{JJ_a}(x) b_{JJ_a}(y) dx dy
%     phi_self = Kxy*(S.W.*gauss);
%     integral = sum(S.W.*phi_self.*gauss);
% 
%     S.Eself = S.Eself + 0.5*integral;
% end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

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

function K = yukawa_kernel(x,y,S)
% Infinite periodic case using cosh and sinh

% Takes two vectors x and y and returns kernel K(x,y)
% Distance matrix with Minimum Image Convention
x = x(:);
y = y(:);
raw_dist = abs(x - y');
dist_matrix = min(raw_dist, S.L - raw_dist);

% Periodic Kernel Formula (Cosh)
% This sums infinite periodic images exactly
numerator = cosh(S.kappa * (0.5 * S.L - dist_matrix));
denominator = 2 * S.kappa * sinh(S.kappa * S.L / 2);

prefactor = 4 * pi / S.epsilon_elec; % Note prefactor change for this form
K = prefactor * (numerator / denominator);
end