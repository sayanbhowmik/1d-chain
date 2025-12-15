function [atompos_next,AtmForce,S] = electronicGroundStateAtomicForce(atom_pos,S)

tic_relax = tic;
fprintf('\n');
fprintf(' ###############################################################\n');
fprintf(' Relaxation step number: %d \n', S.Relax_iter);

% Check position of atom near the boundary and apply wraparound in case of PBC
S = check_atomlocation(S);

% Pseudocharge
t_calc_b = tic;

S = calculate_b_guessRho_eself(S);

% % Self energy
% S = calculate_Eself(S);

fprintf(' Time for b calculation: %.3f seconds.\n',toc(t_calc_b));

% set up guess electron density (guess rho)
S = initElectrondensity(S);

% Self-consistent Field (SCF) method
S = scf(S);

fprintf('\n');
fprintf(' **********************************************************\n');
fprintf(' *          Energy per unit cell = %11.9f Ha.       *\n', S.Etotal);
fprintf(' *          Energy per atom = %11.9f Ha.            *\n', S.Etotal / S.n_atm);
fprintf(' **********************************************************\n');

% Atomic force calculation
tic_force = tic;
S.force = atomicForce(S);
force_mat = S.force;
sz_fmat = size(force_mat);
force_corr = sum(force_mat, 1) / sz_fmat(1);
S.force = force_mat - repmat(force_corr, sz_fmat(1),1);

fprintf(' ***********************************************************\n');
fprintf(' *                      Atomic Force                       *\n');
fprintf(' ***********************************************************\n');
fprintf(' Drift free forces (Ha/Bohr):\n');
disp(S.force);
S.abs_force = sqrt(sum(abs(S.force).^2,2));
fprintf(' Max magnitude of forces (Ha/Bohr):');
disp(max(S.abs_force));

t_force = toc(tic_force);
fprintf('\n Time for calculating forces: %f s.\n', t_force);
fprintf(' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n');
fprintf(' SCF step time                      :  %.3f (sec)\n',toc(tic_relax));
fprintf(' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ \n');

% Assign output of function
AtmForce = reshape(S.force',[],1);
AtmForce = - AtmForce;

% Dummy return for now
atompos_next = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = initElectrondensity(S)
S.rho = S.rho_at;
% perform charge extrapolation : Write later
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function S = check_atomlocation(S)
% map atom positions back to the domain if atoms are outside the domain in
% periodic directions
S.Atoms(S.Atoms(:,1) >= S.L,1) = S.Atoms(S.Atoms(:,1) >= S.L,1) - S.L;
S.Atoms(S.Atoms(:,1) < 0,1) = S.Atoms(S.Atoms(:,1) < 0,1) + S.L;
end