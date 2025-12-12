function S = calculate_b_guessRho(S)
t1 = tic;
fprintf('\n Starting pseudocharge generation...\n');

% initialization
S.b = zeros(S.N,1);
S.Eself = 0;

% Pseudopotential (local)
VJ_mat = zeros(S.N,S.n_typ);
for i = 1 : S.n_typ
    VJ_mat(:,i) = calculate_VJ(S,S.unique_Z(i),S.unique_sigma(i));
end

%--------------------------------------------------------------------------
% We are currently using pseudocharge radii implementation - can be
% switched off if needed later
%--------------------------------------------------------------------------
% Pseudocharge generation
for JJ_a = 1:S.n_atm
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
        S.b(ii_s:ii_e) = S.b(ii_s:ii_e) + bJ(II);
        
        %------------------------------------------------------------------
        %- Numerical implementation -- needs very fine mesh for good
        %- accuracy
        % idx_type = find(S.unique_Z == S.Z(JJ_a));
        % V_PS = interp1(S.x,S.x.*VJ_mat(:,idx_type),abs(dd),'spline');
        % 
        % V_PS = V_PS./abs(dd);
        % V_PS(abs(dd)<S.x(2)) = VJ_mat(1,idx_type);
        % 
        % S.Eself = S.Eself + 0.5*sum(S.dx*(V_PS(II).*bJ(II)));
        %------------------------------------------------------------------

    end % end of loop over all images (including atom JJ_a)

end % end of loop over all atoms

% fprintf('Eself: %.15f \n',S.Eself);

%--------------------------------------------------------------------------
% Naive implementation considering full space
%--------------------------------------------------------------------------

% atom_positions = S.Atoms;
% 
% for i = 1:S.n_atm
%     dist = abs(S.x - atom_positions(i));
%     % Apply Periodic Boundary Condition to distance
%     dist = min(dist, S.L - dist); 
% 
%     gauss = -(S.Z(i) / sqrt(2*pi*S.b_sigma(i)^2)) * exp(-dist.^2 / (2*S.b_sigma(i)^2));
%     S.b = S.b + gauss;
% end
%--------------------------------------------------------------------------

% Verify Charge Neutrality of Background
fprintf(' Total Nuclear Charge (Integral b): %.14f (Target: %.14f)\n', ...
    dot(S.W,S.b), -S.Nelectron);

% find positive charges and negative charges
S.PosCharge = abs(dot(S.W, S.b));
S.NegCharge = -S.PosCharge + S.NetCharge;

% Atomic density guess
S.rho_at = -S.b;

fprintf(' Done. (%f s)\n',toc(t1));
end

%------------------------------------------------------------------
% Useful functions in case of numerical implementation
%------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function VJ = calculate_VJ(S, Z, b_sigma)
    % Create a "Virtual" Extended Grid
    % Extend by 3x to ensure Yukawa tail decays to < 1e-4
    L_virtual = 3 * S.L;
    N_virtual = 3 * S.N; 

    % We need to create a temporary diff_mat for this larger grid
    S_virtual = S;
    S_virtual.L = L_virtual;
    S_virtual.N = N_virtual;
    S_virtual.x = (0:N_virtual-1)' * S.dx; % New long x-axis

    % 2. Setup Matrices for Virtual Grid
    % (Re-use your existing topology generation functions)
    diff_mat = struct('L', L_virtual, 'N', N_virtual, 'dx', S.dx, ...
        'FDn', S.FDn);

    % Finite difference weights of the second derivative
    FDn = diff_mat.FDn;
    w2 = zeros(1,FDn+1) ;
    for k=1:FDn
    	w2(k+1) = (2*(-1)^(k+1))*(factorial(FDn)^2)/...
    		(k*k*factorial(FDn-k)*factorial(FDn+k));
    	w2(1) = w2(1)-2*(1/(k*k));
    end
    diff_mat.w2 = w2;

    % Use your symmetry-corrected Laplacian generator
    diff_mat = lapIndicesValues_1d(diff_mat); 
    DL11_virt = Laplacian_1d(diff_mat);

    x_virt = S_virtual.x;

    % 3. Source (Gaussian on Virtual Grid)
    gauss = -(Z / sqrt(2*pi*b_sigma^2)) * exp(-x_virt.^2 / (2*b_sigma^2));
    RHS = (4*pi/S.epsilon_elec) * gauss;

    % 4. Apply Zero Dirichlet BCs at the FAR end (x = 3L)
    % At 3L, exp(-kappa*3L) ~ exp(-9.6) ~ 6e-5 (Safe to force to 0)
    RHS_reduced = RHS(1:end-1); 

    LHS = -DL11_virt + (S.kappa^2 * speye(N_virtual));
    LHS_reduced = LHS(1:end-1, 1:end-1);

    % 5. Solve
    VJ_virt = LHS_reduced \ RHS_reduced;
    VJ_virt = [VJ_virt; 0]; % Pad back

    % 6. Crop Result back to original size (S.N)
    % We only care about the potential inside the actual box [0, L]
    VJ = VJ_virt(1:S.N);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DL11 = Laplacian_1d(S)
Nx = S.N;

% D_xx laplacian in 1D
%-----------------------
V = S.V_11;

% Create discretized Laplacian
DL11 = sparse(S.I_11,S.II_11,V,Nx,Nx);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = lapIndicesValues_1d(S)
Nx = S.N;
n0 = S.FDn;
w2 = S.w2;
dx = S.dx;

% Initial number of non-zeros: including ghost nodes
nnzCount = (2 * n0 + 1) * Nx;

% Row and column indices and the corresponding non-zero values
% used to generate sparse matrix DL11 s.t. DL11(I(k),II(k)) = V(k)
I = zeros(nnzCount,1);
V = zeros(nnzCount,1);
II = zeros(nnzCount,1);
rowCount = 1;
count = 1;
coef_dxx = 1/dx^2;

% Find non-zero entries that use forward difference
for ii = 1:Nx
    % diagonal element
    I(count) = rowCount; II(count) = ii;
    V(count) = w2(1)*coef_dxx ;
    count = count + 1;
    % off-diagonal elements
    % 2. Off-diagonal elements
    for q = 1:n0
        % --- Right Neighbor (ii + q) ---
        I(count) = rowCount;
        II(count) = ii + q;
        V(count) = w2(q+1) * coef_dxx;
        count = count + 1;

        % --- Left Neighbor (ii - q) ---
        idx_left = ii - q;

        % MIRROR SYMMETRY FIX:
        % If index falls to the left of 1 (indices 0, -1, -2...),
        % reflect it back into the domain.
        if idx_left < 1
            idx_left = 2 - idx_left;
            % Example Check:
            % if idx_left is 0  -> 2 - 0  = 2  (Correct)
            % if idx_left is -1 -> 2 - -1 = 3  (Correct)
        end

        I(count) = rowCount;
        II(count) = idx_left;
        V(count) = w2(q+1) * coef_dxx;
        count = count + 1;
    end
    rowCount = rowCount + 1;
end

% Removing outside domain entries
isIn = (II >= 1) & (II <= Nx);
S.I_11 = I(isIn); S.II_11 = II(isIn); S.V_11 = V(isIn);
end