function [S,x_kp1] = mixing(S, g_k, x_k, iter)
[S,x_kp1] = Periodic_Pulay(S, g_k, x_k, iter);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [S,x_kp1] = Periodic_Pulay(S, g_k, x_k, iter)

m = S.MixingHistory; 
p = S.PulayFrequency;
beta = S.MixingParameter; % beta
omega = S.MixingParameterSimple; % mixing parameter for simple mixing step

Pulay_mixing_flag = (rem(iter,p) == 0 && iter > 1);

if Pulay_mixing_flag   % paulay mixing
    amix = beta; 
else                   % simple (linear) mixing
    amix = omega; 
end

f_k = g_k - x_k;
if iter > 1
	f_km1 = S.mixing_hist_fkm1;
	x_km1 = S.mixing_hist_xkm1;
end

% store residual & iteration history
if iter > 1
	i_hist = mod(iter-2,m)+1;
	if (S.PulayRestartFlag ~= 0 && i_hist == 1)
		S.X = zeros(size(S.X)); S.F = zeros(size(S.F));
		S.X(:,1) = x_k - x_km1;
		S.F(:,1) = f_k - f_km1;
	else
		S.X(:,i_hist) = x_k - x_km1;
		S.F(:,i_hist) = f_k - f_km1;
	end
end

% apply Anderson extrapolation every p iters
if Pulay_mixing_flag
    % find weighted averages x_wavg, f_wavg, Gamma
    [x_wavg, f_wavg, Gamma] = andersonWtdAvg(x_k(1:S.nspden*S.N), f_k(1:S.nspden*S.N), ...
        S.X(1:S.nspden*S.N,:), S.F(1:S.nspden*S.N,:),S.nspden,S.MixingVariable);
    S.mix_Gamma = Gamma; % store Gamma to use in Hubbard DFT if required
else 
	% simple mixing
	x_wavg = x_k; f_wavg = f_k;
end

% calculate sum of all columns
f_wavg = reshape(f_wavg,[],S.nspden);

Pf = zeros(S.N,S.nspden);
Pf(:,1) = amix * f_wavg(:,1);     % mixing param is included in Pf, no preconditioner

% Flatten Pf
Pf = reshape(Pf,[],1);

% get x_kp1
x_kp1 = x_wavg + reshape(Pf,[],1);

negrho_count = sum(x_kp1(1:S.N) < 0);
if (negrho_count > 0)
    fprintf('\nDensity got negative\n\n');                
end
x_kp1(x_kp1(1:S.N) < 0) = 0;
integral = S.W'*(x_kp1(1:S.N));
x_kp1(1:S.N) = x_kp1(1:S.N) * (-S.NegCharge/integral);

% update the history vectors
S.mixing_hist_fkm1 = f_k;
S.mixing_hist_xkm1 = x_k;
S.mixing_hist_xk = x_kp1;

end