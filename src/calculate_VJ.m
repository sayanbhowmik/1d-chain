%------------------------------------------------------------------
% Analytical VJ to bJ
%------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function VJ = calculate_VJ(x, S, JJ_a)
    % Inputs:
    % x     : Coordinate vector
    % S     : Struct containing .kappa and .epsilon_elec
    % JJ_a  : Index of atom
    
    % Atomic number and sigma
    Z = S.unique_Z(JJ_a); sigma = S.unique_sigma(JJ_a);

    % 1. Constants
    k = S.kappa;
    eps = S.epsilon_elec;
    
    % The prefactor derived: - (pi * Z) / (eps * k)
    % Note: The '4' from 4pi canceled with the denominators during integration.
    prefactor = - (pi * Z) / (eps * k);
    
    % The exponential scaling factor
    exp_factor = exp((k^2 * sigma^2) / 2);
    
    % 2. Error Function Arguments
    sqrt2_sig = sqrt(2) * sigma;
    
    % Argument 1: (k*sigma^2 - x) / (sqrt(2)*sigma)
    arg1 = (k * sigma^2 - x) / sqrt2_sig;
    
    % Argument 2: (k*sigma^2 + x) / (sqrt(2)*sigma)
    arg2 = (k * sigma^2 + x) / sqrt2_sig;
    
    % 3. Calculate Terms
    % term1 corresponds to the integral from -inf to x
    term1 = exp(-k * x) .* erfc(arg1);
    
    % term2 corresponds to the integral from x to inf
    term2 = exp( k * x) .* erfc(arg2);
    
    % 4. Combine
    VJ = prefactor * exp_factor * (term1 + term2);
end