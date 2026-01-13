function S = exx_initialization(S)
% ACE FLAG - Default on
S.ACEFlag = 1; % 1 - Use ACE, 0 - Do not use ACE (not optimized at all for now)

% Method - Fourier vs Real
S.exxmethod = 0; % 0 - Fourier, 1 - Real

S.MAXIT_FOCK = 30;
S.MINIT_FOCK = 2;
S.FOCK_TOL = 0.2*S.SCF_tol;
S.SCF_tol_init = max(10*S.FOCK_TOL,1e-3);
S.hyb_mixing = 0.25;
S.hyb_range_fock = 0.1587;          % VASP
S.hyb_range_pbe = 0.1587;           % VASP
S.EXXACEVal_state = 3;

S = const_for_FFT(S);

if S.exxmethod == 1 || S.ACEFlag == 0
    % Define Yukawa LHS
    LHS_yuk = -S.Lap_std + (S.kappa^2 * speye(S.N));

    % Pre-compute Screening Operator (For real space)
    exponent_arg = -0.25 / (S.hyb_range_fock^2) * full(LHS_yuk);
    S.Screening_Op = eye(S.N) - expm(exponent_arg);

    % Pre-factorize the Solver for speed
    % decomposition object caches the Cholesky factor
    S.dA = decomposition(LHS_yuk, 'chol');

end

end

function S = const_for_FFT(S)
    % Setup Dimensions
    N = S.N;
    L = N * S.dx;
    
    % Parameters
    mu = S.hyb_range_fock;  % HSE Screening (Omega)
    kappa2 = S.kappa^2;     % Yukawa Screening (from Image)
    
    % Dielectric Constant
    if isfield(S, 'epsilon_elec')
        eps_val = S.epsilon_elec;
    else
        eps_val = 1.0; 
        warning('S.epsilon_elec not found. Defaulting to 1.0');
    end

    S.const_by_alpha = zeros(1, N);
    S.k_shift = [0, 0, 0];
    S.num_shift = 1;

    % Vectorized G^2
    k_indices = [0:floor(N/2), -floor(N/2)+1-mod(N,2):-1];
    G_squared = (k_indices * 2 * pi / L).^2;
    
    % The EXACT Hybrid Kernel
    % Denominator: Matches the Yukawa Green's Function (-d2/dx2 + k2)
    denom = G_squared + kappa2;
    
    % Exponent: Matches the Gaussian Screening (erfc)
    % Note: We screen based on the TOTAL energy scale (G^2 + k^2)
    arg = denom / (4 * mu^2);
    
    % Combine them
    % The (1 - exp) term handles the short-range HSE cut-off
    kernel = (4 * pi ./ (eps_val * denom)) .* (1 - exp(-arg));
    
    % Zero out G=0 for stability 
    % This removes the "background charge" energy shift (needed if kappa = 0)
    if sqrt(kappa2) == 0
        kernel(1) = 0;
    end
    
    S.const_by_alpha(1, :) = kernel;
end

