function [S] = scf(S)
fprintf('\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n');
fprintf(' Starting SCF iteration...\n');

% Electrostatic potential
S = poissonSolve(S);

% Exchange-correlation potential
S = exchangeCorrelationPotential(S);

% Effective potential
S = calculate_effective_potential(S);

% Initialize the mixing history vectors
S.X = zeros(S.N*S.nspden,S.MixingHistory);
S.F = zeros(S.N*S.nspden,S.MixingHistory);
S.mixing_hist_xkm1 = zeros(S.N*S.nspden,1);
S.mixing_hist_fkm1 = zeros(S.N*S.nspden,1);

% initialize history
S.mixing_hist_xkm1(1:S.N) = S.rho(:,1);
S.mixing_hist_xk = S.mixing_hist_xkm1;

S.lambda_f = 0.0;
S = scf_loop(S,S.SCF_tol);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = scf_loop(varargin)
S = varargin{1};
scf_tol = varargin{2};

% SCF LOOP 
count_SCF = 1;

% time for one scf
t_SCF = 0; 

% start scf loop
while count_SCF <= S.MAXIT_SCF
    tic_eig = tic;

    rho_in = S.rho; % For mixing

    fprintf(' ========================= \n');
    fprintf(' Relaxation iteration: %2d \n SCF iteration number: %2d \n',S.Relax_iter,count_SCF);
    fprintf(' ========================= \n');

    [S.psi, S.EigVal] = eigSolver(S);

    % Solve for Fermi energy S.lambda_f and occupations
    S = occupations(S);

    % for density mixing, us rho_in to estimate total energy
    [S.Etotal,S.Eband,S.Exc,S.Exc_dc,S.Eelec_dc,S.Eent] = evaluateTotalEnergy(S);

    % Electron density
    S = electronDensity(S);

    % Error in SCF fixed-point iteration
    err = (norm(S.rho - rho_in))/(norm(S.rho));
    fprintf(' Error in SCF iteration: %.4e \n',err) ;

    % Mixing to accelerate SCF convergence
    % get rho_out
    rho_out = zeros(S.N*S.nspden,1);
    rho_out(1:S.N) = S.rho(:,1);
    % mixing
    [S, rho_out] = mixing(S,rho_out,S.mixing_hist_xk,count_SCF);
    % update
    S.rho(:,1) = rho_out(1:S.N);

    % update Veff
    S = poissonSolve(S);
    % Exchange-correlation potential
    S = exchangeCorrelationPotential(S);
    % Effective potential            
    S = calculate_effective_potential(S);

    scf_runtime = toc(tic_eig);
	t_SCF = t_SCF + scf_runtime; % add chebyshev filtering time to SCF time
    fprintf(' This SCF iteration took %.3f s.\n\n', scf_runtime);

    count_SCF = count_SCF + 1 ;

    if err < scf_tol && count_SCF > S.MINIT_SCF
        break;
    end
end

fprintf('\n Finished SCF iteration in %d steps!\n', (count_SCF - 1));
fprintf(' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [S] = calculate_effective_potential(S)
S.Veff = S.XCswitch*S.Vxc + S.phi;
S.Veff = real(S.Veff);
end