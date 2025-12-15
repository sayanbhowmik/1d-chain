function S = calculate_b_guessRho_eself(S)
t1 = tic;
fprintf('\n Starting pseudocharge generation...\n');

% initialization
S.b = zeros(S.N,1);
S.Eself = 0;

% Pseudopotential (local)
VJ_mat = zeros(S.N,S.n_typ);
for i = 1 : S.n_typ
    VJ_mat(:,i) = calculate_VJ(S.x,S,i);
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

        % Pseudopotential at grid points through interpolation
        idx_type = find(S.unique_Z == S.Z(JJ_a));
        V_PS = interp1(S.x,VJ_mat(:,idx_type),abs(dd),'spline');

        % Pseudocharge density calculation - numerical bJ
        II = 1+S.FDn : size(dd,1)-S.FDn;
        bJ = pseudochargeDensity_atom(V_PS,II,S);
        S.b(ii_s:ii_e) = S.b(ii_s:ii_e) + bJ(II);

        % Eself
        S.Eself = S.Eself + 0.5*sum(S.dx*(V_PS(II).*bJ(II)));

        % % Pseudocharge density calculation - analytical bJ
		% II = 1+S.FDn : size(dd,1)-S.FDn;
        % 
        % bJ = pseudochargeDensity_atom(dd,II,S.Z(JJ_a),S.b_sigma(JJ_a));
        % S.b(ii_s:ii_e) = S.b(ii_s:ii_e) + bJ(II);
        
        %------------------------------------------------------------------
        %- Numerical implementation Eself -- needs very fine mesh for good
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

fprintf(' Eself: %.15f \n',S.Eself);

% %--------------------------------------------------------------------------
% % Naive implementation considering full space
% %--------------------------------------------------------------------------
% 
% atom_positions = S.Atoms;
% b = zeros(S.N,1);
% for i = 1:S.n_atm
%     dist = abs(S.x - atom_positions(i));
%     % Apply Periodic Boundary Condition to distance
%     dist = min(dist, S.L - dist); 
% 
%     gauss = -(S.Z(i) / sqrt(2*pi*S.b_sigma(i)^2)) * exp(-dist.^2 / (2*S.b_sigma(i)^2));
%     % S.b = S.b + gauss;
%     b = b + gauss;
% end
% fprintf(' Normed difference ||b_num - b_analytic|| = %.14f\n',norm(S.b-b));
% %--------------------------------------------------------------------------

% Verify Charge Neutrality of Background
fprintf(' Total Nuclear Charge (Integral b): %.14f (Target: %.14f)\n', ...
    dot(S.W,S.b), -S.Nelectron);

% find positive charges and negative charges
S.PosCharge = abs(dot(S.W, S.b));
S.NegCharge = -S.PosCharge + S.NetCharge;

% % Atomic density guess - different
% atom_positions = S.Atoms;
% S.rho_at = zeros(S.N,1);
% for i = 1:S.n_atm
%     dist = abs(S.x - atom_positions(i));
%     % Apply Periodic Boundary Condition to distance
%     dist = min(dist, S.L - dist); 
%     sigma = 1.0*S.b_sigma(i);
% 
%     gauss = (S.Z(i) / sqrt(2*pi*sigma^2)) * exp(-dist.^2 / (2*sigma^2));
%     % S.b = S.b + gauss;
%     S.rho_at = S.rho_at + gauss;
% end
% rho_scal = abs(S.NegCharge/dot(S.W,S.rho_at));
% S.rho_at = rho_scal*S.rho_at;
% fprintf(' Integral rho = %0.14f\n',dot(S.W,S.rho_at));

% Atomic density guess
S.rho_at = -S.b;

fprintf(' Done. (%f s)\n',toc(t1));
end