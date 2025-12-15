%--------------------------------------------------------------------------
% Numerical b
%--------------------------------------------------------------------------
function b = pseudochargeDensity_atom(V,II,S)
% Calculate [epsilon/(4pi)] [-lap * V + k^2 * V] = b (in this function)
b = zeros(size(V));

% lap*V
dx2 = S.dx*S.dx;
coeff =  S.w2(1) * (1/dx2);
b(II) = coeff * V(II);
for p = 1:S.FDn
    b(II) = b(II) + S.w2(p+1)/dx2 * (V(II+p) + V(II-p));
end
b = -b; % -lap*V

% k^2*V
coeff = S.kappa^2;
b(II) = b(II) + coeff * V(II);

% Multiply with (epsilon/(4pi))
b(II) = (S.epsilon_elec/(4*pi))*b(II);
end

%--------------------------------------------------------------------------
% Analytic b
%--------------------------------------------------------------------------
% function b = pseudochargeDensity_atom(dd,II,Z,sigma)
%
% % % For now only one type of atom
% % sigma = S.b_sigma(idx);
% % Z = S.Z(idx);
%
% b = zeros(size(dd));
%
% b(II) = (-Z/sqrt(2*pi*sigma^2))*exp(-(0.5/(sigma^2)) * (dd(II).^2) );
% end