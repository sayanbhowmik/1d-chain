function force = atomicForce(S)

% Force_J = -dEtotal/dR_J
% Force_J = -dE_el/dR_J + dE_self/dR_J

fprintf('\n Starting atomic force calculation ... \n');

% Calculate local forces
tic_locforces = tic;

force_local = zeros(S.n_atm,1);
force_self = zeros(S.n_atm,1);

% % Calculate Yukawa Kernel
% Kxy = yukawa_kernel(S.x,S.x,S);

% Pseudopotential (local)
VJ_mat = zeros(S.N,S.n_typ);
for i = 1 : S.n_typ
    VJ_mat(:,i) = calculate_VJ(S.x,S,i);
end

Dphi_x = S.grad*(S.phi); % gradient of phi

unique_rows = [S.unique_Z, S.unique_sigma];
for JJ_a = 1 : S.n_atm
    % Initialize b_j and db_j/dx on full grid size
    % bJ_full = zeros(S.N,1);
    dbJ_full = zeros(S.N,1);

    % Atom position of atom JJ_a
    x0 = S.Atoms(JJ_a,1);

    % Find index of rb_x for the atom type
    % idx_type = find(S.unique_Z == S.Z(JJ_a));
    row = [S.Z(JJ_a), S.b_sigma(JJ_a)];
    idx_type = find( all(unique_rows == row, 2));

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

        % Pseudopotential at grid points through interpolation
        % idx_type = find(S.unique_Z == S.Z(JJ_a));
		V_PS = interp1(S.x,VJ_mat(:,idx_type),abs(dd),'spline');

        % Pseudocharge density calculation - numerical bJ
        II = 1+S.FDn : size(dd,1)-S.FDn;
        bJ = pseudochargeDensity_atom(V_PS,II,S);

        dVJ_x = dpseudopot(V_PS,II,S);

        II_rb = ii_s : ii_e;

        fl_x = -(sum(S.W(II_rb).*bJ(II).*(Dphi_x(II_rb) - dVJ_x(II))));
        force_local(JJ_a,1) = force_local(JJ_a,1) + fl_x;
        force_self(JJ_a,1) = force_self(JJ_a,1) -(sum(S.W(II_rb).*bJ(II).*(-dVJ_x(II))));

        % % Pseudocharge density calculation - analytical
		% II = 1+S.FDn : size(dd,1)-S.FDn;
        % 
        % bJ = pseudochargeDensity_atom(dd,II,S.Z(JJ_a),S.b_sigma(JJ_a));
        % dbJ = -(dd/S.b_sigma(JJ_a)^2).*bJ;
        % 
        % II_act = ii_s:ii_e;
        % 
        % % bJ_full(II_act) = bJ_full(II_act) + bJ(II);
        % dbJ_full(II_act) = dbJ_full(II_act) + dbJ(II);
        % 
        % phi_bj = analytic_yukawa_solution(dd, S, JJ_a);
        % fl_x2 = -sum(S.W(II_act).*dbJ(II).*phi_bj(II));
        % 
        % force_local(JJ_a,1) = force_local(JJ_a,1) + fl_x2;
        % force_self(JJ_a,1) = force_self(JJ_a,1) + fl_x2;
    end
    % fl_x = sum(S.W.*dbJ_full.*S.phi);
    % force_local(JJ_a,1) = force_local(JJ_a,1) + fl_x;

    % phi_self = Kxy*(S.W.*bJ_full); % Yukawa kernel way

    % % Poisson way
    % RHS = (4*pi/S.epsilon_elec)*(bJ_full);
    % A = -S.Lap_std + (S.kappa^2 * speye(S.N));
    % 
    % phi_self = A \ RHS;
    % 
    % fl_x2 = -sum(S.W.*dbJ_full.*phi_self);
    % 
    % force_local(JJ_a,1) = force_local(JJ_a,1) + fl_x2;
    % force_self(JJ_a,1) = force_self(JJ_a,1) + fl_x2;
end
fprintf(' local force calculation: %.3f s\n', toc(tic_locforces));

force = force_local;
% fprintf("%0.14f\n",force_self); % test - self force should be 0 if not
% for eggbox

end

function [DX_x] = dpseudopot(X,II,S)
DX_x = zeros(size(X));
for p = 1:S.FDn
    DX_x(II) = DX_x(II) + S.w1(p+1)/S.dx*(X(II+p)-X(II-p));
end
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