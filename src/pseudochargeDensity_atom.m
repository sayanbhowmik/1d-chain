function b = pseudochargeDensity_atom(dd,II,Z,sigma)

% % For now only one type of atom
% sigma = S.b_sigma(idx);
% Z = S.Z(idx);

b = zeros(size(dd));

b(II) = (-Z/sqrt(2*pi*sigma^2))*exp(-(0.5/(sigma^2)) * (dd(II).^2) );
end