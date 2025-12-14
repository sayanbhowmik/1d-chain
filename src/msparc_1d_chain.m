function S = msparc_1d_chain(atm_dist, n_atm, kappa, epsilon, XC)
% function S = msparc_1d_chain()
if strcmpi(XC,'None')
    XCswitch = 0; % Enable exchange-correlation functional
    XC = 'GGA_PBE';
else
    XCswitch = 1; % Disable exchange-correlation functional
end
format long;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: Consider all units to be atomic units unless otherwise stated
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start timer
total_time = tic;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%          Atomic information             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Equidistant atoms
% atm_dist = 10; % interatomic distance, user specified
% S.n_atm = floor(S.L/atm_dist);
Atoms = (0:n_atm-1)' *atm_dist;
% S.Atoms = (0:S.n_atm-1)' * atm_dist; % Column vector of positions [0, a, 2a, ...]

% Non-Equidistant atoms - define as you like
% atm_dist = 10; % interatomic distance, user specified
% n_atm = 28;
% S.n_atm = 16;
% S.n_atm = n_atm;
% Atoms = (0:n_atm-1)' *atm_dist;
% S.Atoms = (0:S.n_atm-1)' * atm_dist; % Column vector of positions [0, a, 2a, ...]

% % Force numerical test - perturb any atom
% S.Atoms(end) = S.Atoms(end) - 0.1;
% Atoms(end) = Atoms(end) - 0.1;

% Store atomic numbers and b_sigma as a vector - change accordingly
% S.Z = 2*ones(size(S.Atoms)); 
Z = 2*ones(size(Atoms));
%S.b_sigma = 2*ones(size(S.Atoms)); % sigma of pseudocharge gaussian in Lin Lin paper
b_sigma = 2*ones(size(Atoms)); % sigma of pseudocharge gaussian

% S.Nelectron = sum(S.Z);
Nelectron = sum(Z); % Total number of electrons

% Find unique pairs of (Z, sigma)
% [unique_rows, ~, ~] = unique([S.Z(:), S.b_sigma(:)], 'rows', 'stable');
unique_rows = unique([Z(:), b_sigma(:)], 'rows', 'stable');
unique_Z = unique_rows(:, 1);
unique_sigma = unique_rows(:, 2);

% S.n_typ = length(unique_Z);
% S.unique_Z = unique_Z;
% S.unique_sigma = unique_sigma;
n_typ = length(unique_Z);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basic parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L = 320; % Lattice size
L = atm_dist * n_atm
% N = 3200; % Number of grid points
N = floor(10 * L); % Number of grid points based on lattice size
dx = L/N; % grid spacing
SCF_tol = 1e-6; % SCF tolerance

% FDn: half of finite difference order
FDn = 6;

S = struct('L',L,'N',N,'dx',dx,'SCF_tol',SCF_tol,'FDn',FDn,'n_atm',n_atm,'Atoms',Atoms,'Z',Z,'b_sigma',b_sigma,'Nelectron',Nelectron,'n_typ',n_typ,'unique_Z',unique_Z,'unique_sigma',unique_sigma);

% Mixing parameter - Anderson mixing only (no preconditioners)
% Setting default to 0.3, reduce to 0.1 in case of difficulty
S.MixingParameter = 0.1;

% Relax Flag
S.RelaxFlag = 0;
S.Relax_iter = 1;

% XC: Exchange-correlation functional 
% S.XCswitch = 1; % 0 to switch off, 1 to switch on
S.XCswitch = XCswitch;
% S.XC = 'GGA_PBE';
S.XC = XC;
S.isgradient = 0; % default
% decomposition of XC, ixc = [iexch,icorr imeta ivdw]
if strcmp(S.XC, 'LDA_PW')
	S.xc = 0;
    S.ixc = [1 2 0 0];
elseif strcmp(S.XC, 'LDA_PZ')
	S.xc = 1; 
    S.ixc = [1 1 0 0];
elseif strcmp(S.XC, 'GGA_PBE')
	S.xc = 2;
    S.ixc = [2 3 0 0];
    S.xc_option = [1 1];
    S.isgradient = 1;
elseif strcmp(S.XC, 'GGA_PBEsol')
    S.ixc = [2 3 0 0];
    S.xc_option = [2 2];
    S.isgradient = 1;
elseif strcmp(S.XC, 'GGA_RPBE')
    S.ixc = [2 3 0 0];
    S.xc_option = [3 3];
    S.isgradient = 1;
end

% Yukawa parameter for electrostatics
% S.kappa = 0.01;
S.kappa = kappa;
% S.epsilon_elec = 10;
S.epsilon_elec = epsilon;

% % Electronic smearing - specify width
% S.elec_T_type = 0; % fermi-dirac
% smearing = 0.1; % eV
% S.bet = 27.21138602 / smearing; % smearing = 0.1 eV = 0.00367493225 Ha, Beta := 1 / smearing
% S.Temp = 1./(3.166810501187400e-06 * S.bet);

% Electronic smearing - specify Temp
% Cst: Factor for conversion from Ha to eV
Cst = 27.21138602;
S.kB = (8.6173303e-5)/Cst;
S.elec_T_type = 0; % fermi-dirac
S.Temp = 100; % 100 K
S.bet = 1 / (S.kB*S.Temp);

% Initialize remaining parameters
S = initialization(S);

% Perform relaxation/single point calculation
if S.RelaxFlag 
    error(" Relaxation not implemented.\n");
else
    [~,~,S] = electronicGroundStateAtomicForce(S.Atoms,S);
end

t_wall = toc(total_time);
% Program run-time
fprintf('\n Run-time of the program: %f seconds\n', t_wall);
end