function S = initialization(S)
S.x = (0:S.N-1)' * S.dx;  % Grid vector

% SCF
S.NetCharge = 0;
S.MAXIT_SCF = 100;
S.MINIT_SCF = 1;
S.TOL_LANCZOS = 0.01*S.SCF_tol;

% Density tolerance for exchange-correlation
S.xc_rhotol = 1e-14;
S.xc_magtol = 1e-8;
S.xc_sigmatol = 1e-24;

% Number of spin density channels
S.nspden = 1;
S.nspin = 1;
S.nspinor = 1;
% S.occfac = 2/S.nspinor;
S.occfac = 1;
S.nspinor_eig = S.nspinor/S.nspin;

% Nev - Number of states
fprintf('## Finding Nev ...\n');
% S.Nev = S.nspinor_eig*(floor(S.Nelectron / 2) * 1.2 + 5); 
S.Nev = S.nspinor_eig*(floor(S.Nelectron / 2) * 2 + 10); 
S.Nev = round(S.Nev);
fprintf('## Based on the number of electrons, Nev is set to: %d\n',S.Nev);

% Mixing
S.MixingVariable = 0;
S.MixingHistory = 7;
S.MixingParameterSimple = S.MixingParameter;

S.PulayFrequency = 1;
S.PulayRestartFlag = 0;

% Preconditiong - Kerker
S.MixingPrecond = 1; % Default with Kerker preconditioning
S.precond_kerker_kTF = 1;
S.precond_kerker_thresh = 0.1;

% Weights for spatial integration over domain
S.W = ones(S.N,1)*S.dx;

% Finite difference weights of the second derivative
FDn = S.FDn;
w2 = zeros(1,FDn+1) ;
for k=1:FDn
	w2(k+1) = (2*(-1)^(k+1))*(factorial(FDn)^2)/...
		(k*k*factorial(FDn-k)*factorial(FDn+k));
	w2(1) = w2(1)-2*(1/(k*k));
end
S.w2 = w2;

fprintf(' Creating differentiation matrices ...\n');
% Calculate discrete laplacian for periodic only
S = lapIndicesValues_1d(S);
S.Lap_std = Laplacian_1d(S);
% S.Lap_std = (S.Lap_std + S.Lap_std') / 2; % Make it symmetric

% Calculate preconditioners for negative discrete laplacian
[S.LapPreconL, S.LapPreconU] = ilu(S.Lap_std,struct('droptol',1e-5));

% Finite difference weights of the first derivative
w1 = zeros(1,FDn) ;
for k=1:FDn
	w1(k+1) = ((-1)^(k+1))*(factorial(FDn)^2)/...
		(k*factorial(FDn-k)*factorial(FDn+k));
end
S.w1 = w1;

% Calculate gradient
S = gradIndicesValues(S);
S.grad = Gradient(S);

%--------------------------------------------------------------------------
% We are currently using pseudocharge radii implementation - can be
% switched off if needed later.
%--------------------------------------------------------------------------
% pseudocharge tol
S.pseudocharge_tol = 0.01*S.SCF_tol;

% Calculate rb
S = Calculate_rb(S);

end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = Calculate_rb(S)
pos_atm_x = 0; % atom location in x-direction
% rb_up_x = (S.dx < 1.5) * (10+10*S.dx) + (S.dx >=1.5) * (20*S.dx-9.5);
rb_up_x = 5*((S.dx < 1.5) * (10+10*S.dx) + (S.dx >=1.5) * (20*S.dx-9.5));

ii_s_temp = -ceil(rb_up_x/S.dx);
ii_e_temp = ceil(rb_up_x/S.dx);

xx_temp = pos_atm_x + (ii_s_temp-S.FDn:ii_e_temp+S.FDn)*S.dx;
Nx = (ii_e_temp-ii_s_temp)+1;

dd_temp = bsxfun(@minus,xx_temp,pos_atm_x)';
% dd_temp = sqrt(dd_temp.^2);

% Dirichlet integration weights
W_temp = ones(Nx,1)*S.dx;
W_temp(1) = W_temp(1)*0.5; W_temp(Nx) = W_temp(Nx)*0.5;

% Pseudopotential (local)
VJ_mat = zeros(S.N,S.n_typ);
for i = 1 : S.n_typ
    VJ_mat(:,i) = calculate_VJ(S.x,S,i);
end

% Initialize rb_x
S.rb_x = zeros(S.n_typ,1);
for ityp = 1:S.n_typ
    II_temp = 1+S.FDn : size(dd_temp,1)-S.FDn;
    
    % % Analytical b
    % b_temp = pseudochargeDensity_atom(dd_temp,II_temp,S.unique_Z(ityp),S.unique_sigma(ityp));

    % Numerical b
    % Pseudopotential at grid points through interpolation
    V_PS = interp1(S.x,VJ_mat(:,ityp),abs(dd_temp),'spline');
    b_temp = pseudochargeDensity_atom(V_PS,II_temp,S);

    rb_x = S.unique_sigma(ityp);
    rb_x = ceil(rb_x/S.dx-1e-12)*S.dx;
    err_rb = 100;
    count = 1;
    fprintf(' Finding rb ...\n');
    while (err_rb > S.pseudocharge_tol && count <= 1000 && rb_x <= rb_up_x)
        rb_x = rb_x + S.dx;
        ii_rb = -1*ii_s_temp+S.FDn-floor(rb_x/S.dx)+1:-1*ii_s_temp+S.FDn+floor(rb_x/S.dx)+1;
        err_rb = abs( sum(W_temp(ii_rb-S.FDn).*b_temp(ii_rb)) + S.unique_Z(ityp) );
        fprintf(' rb = {%0.3f}, int_b = %0.15f, err_rb = %.3e\n',...
            rb_x,sum(W_temp(ii_rb-S.FDn).*b_temp(ii_rb)), err_rb);
        count = count + 1;
    end
    assert(rb_x<=rb_up_x,'Need to increase upper bound for rb!');
    S.rb_x(ityp) = rb_x;
    fprintf('rb = {%.3f} for Type: %d\n',S.rb_x(ityp),ityp);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function DG = Gradient(S)

I = S.G_I_1;
J = S.G_J_1;
V = S.G_V_1;
DG = sparse(I,J,V,S.N,S.N);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = gradIndicesValues(S)
Nx = S.N;
n0 = S.FDn;
w1 = S.w1;
dx = S.dx;

% Initial number of non-zeros: including ghost nodes
nnz_count = 2*n0*Nx ;

% Row numbers and non-zero values
I = zeros(nnz_count,1) ;
V = zeros(nnz_count,1) ;

% Indices of the columns
II = zeros(nnz_count,1);

% Gradient along x_direction
row_count = 1;
count = 1 ;

for ii=1:Nx
	% off-diagonal elements
	for p=1:n0
		% ii+p
		I(count) = row_count; II(count) = ii+p; 
		V(count) = w1(p+1)/dx;
		count = count + 1;
		% ii-p
		I(count) = row_count; II(count) = ii-p ; 
		V(count) = -w1(p+1)/dx;
		count = count + 1;
	end
	row_count = row_count+1;
end

S.G_JOutl_1 = (II<1); S.G_JOutr_1 = (II>Nx); % Warning: Assumed influence of only neighboring cells
II = mod(II+(Nx-1),Nx)+1;

% Getting linear indices of the columns
S.G_J_1 = II;
S.G_I_1 = I;
S.G_V_1 = V;

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
	for q = 1:n0
		% ii + q
		I(count) = rowCount; II(count) = ii+q;
		V(count) = w2(q+1)*coef_dxx;
		count = count + 1;
		% ii - q
		I(count) = rowCount; II(count) = ii-q;
		V(count) = w2(q+1)*coef_dxx;
		count = count + 1;
		
	end
	rowCount = rowCount + 1;
end

S.isOutl_11 = (II<1); S.isOutr_11 = (II>Nx); % Warning: Assumed influence of only neighboring cells
S.I_11 = I; S.II_11 = mod(II+(Nx-1),Nx)+1; S.V_11 = V;
end