function [S] = exchangeCorrelationPotential(S)
rho = S.rho;
rho(rho < S.xc_rhotol) = S.xc_rhotol;
if S.isgradient
    drho = S.grad * rho;
    sigma = drho.*drho;
    sigma(sigma < S.xc_rhotol) = S.xc_rhotol;
end

% iexch 
switch S.ixc(1)
    case 1
        [ex,vx] = slater(rho);
        v2x = zeros(size(rho));
    case 2
        [ex,vx,v2x] = pbex(rho,sigma,S.xc_option(1));
    otherwise
        ex = zeros(size(rho));
        vx = zeros(size(rho));
        v2x = zeros(size(rho));
end

% icorr
switch S.ixc(2)
    case 1
        [ec,vc] = pz(rho);
        v2c = zeros(size(rho));
    case 2
        [ec,vc] = pw(rho);
        v2c = zeros(size(rho));
    case 3
        [ec,vc,v2c] = pbec(rho,sigma,S.xc_option(2));
    otherwise
        ec = zeros(size(rho));
        vc = zeros(size(rho));
        v2c = zeros(size(rho));
end    

exc = ex + ec;
vxc = vx + vc;
v2xc = v2x + v2c;

if S.isgradient % xc involves gradient of rho
    vxc = vxc - S.grad * (v2xc.*drho);
end

S.e_xc = exc;
S.Vxc = vxc;
S.dvxcdgrho = v2xc;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ex,vx] = slater(rho)
% slater exchange
% @param rho  = total electron density

% parameter 
C2 = 0.73855876638202; % 3/4*(3/pi)^(1/3)
C3 = 0.9847450218427;  % (3/pi)^(1/3)

% computation
ex = - C2 * rho.^(1./3.);
vx = - C3 * rho.^(1./3.);
end

function [ec,vc] = pw(rho)
% pw correlation
% @param rho  = total electron density
% J.P. Perdew and Y. Wang, PRB 45, 13244 (1992)

% parameters
A = 0.031091 ;
alpha1 = 0.21370 ;
beta1 = 7.5957 ;
beta2 = 3.5876 ;
beta3 = 1.6382 ;
beta4 = 0.49294 ;

% computation
rs = (0.75./(pi*rho)).^(1./3.);
rsm12 = rs.^(-0.5);
rs12 = rs.^0.5;
rs32 = rs.^1.5;
rs2 = rs.^2;

om = 2*A*(beta1*rs12 + beta2*rs + beta3*rs32 + beta4*rs2);
dom = A*(beta1*rsm12+ 2*beta2 + 3*beta3*rs12 + 2*2*beta4*rs);
olog = log(1 + 1./om);
t = -2*A*(1+alpha1*rs);
ec = t.*olog;
vc = ec - (rs/3.).*(-2*A*alpha1*olog - (t.*dom)./(om.*om+om) ) ;
end

function [ec,vc] = pz(rho)
% pz correlation
% @param rho  = total electron density
% J.P. Perdew and A. Zunger, PRB 23, 5048 (1981).

% parameters
A = 0.0311;
B = -0.048;
C = 0.002;
D = -0.0116;
gamma1 = -0.1423;
beta1 = 1.0529;
beta2 = 0.3334;

% compuatation
ec = zeros(size(rho,1),1);
vc = zeros(size(rho,1),1);
rs = (0.75./(pi*rho)).^(1.0/3.0) ;
islt1 = (rs < 1.0);
lnrs = log(rs(islt1));
sqrtrs = sqrt(rs(~islt1));
ec(islt1) = A * lnrs + B + C * rs(islt1) .* lnrs + D * rs(islt1);
ox = 1.0 + beta1*sqrtrs + beta2*rs(~islt1);
ec(~islt1) = gamma1 ./ ox;
vc(islt1) = lnrs.*(A + (2.0/3.0)*C*rs(islt1)) + (B-(1.0/3.0)*A) + (1.0/3.0)*(2.0*D-C)* rs(islt1);
vc(~islt1) = ec(~islt1) .* (1 + (7.0/6.0)*beta1*sqrtrs + (4.0/3.0)*beta2*rs(~islt1)) ./ ox;
end

function [ex,v1x,v2x] = pbex(rho,sigma,iflag)
% pbe exchange:
% @param rho  = total electron density
% @param grho = |\nabla rho|^2
% @param iflag options
% iflag=1  J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996)
% iflag=2  PBEsol: J.P.Perdew et al., PRL 100, 136406 (2008)
% iflag=3  RPBE: B. Hammer, et al., Phys. Rev. B 59, 7413 (1999)
% iflag=4  Zhang-Yang Revised PBE: Y. Zhang and W. Yang., Phys. Rev. Lett. 80, 890 (1998)
assert(iflag == 1 || iflag == 2 || iflag == 3 || iflag == 4);

% parameters 
mu_ = [0.2195149727645171 10.0/81.0 0.2195149727645171  0.2195149727645171];
mu = mu_(iflag);
kappa_ = [0.804 0.804 0.804 1.245];
kappa = kappa_(iflag);
threefourth_divpi = 3.0/4.0/pi;
sixpi2_1_3 = (6.0 * pi^2)^(1.0/3.0);
sixpi2m1_3 = 1.0/sixpi2_1_3;
mu_divkappa = mu/kappa;

% computation
rho_updn = rho/2.0;
rho_updnm1_3 = rho_updn.^(-1.0/3.0);
rhomot = rho_updnm1_3;
ex_lsd = -threefourth_divpi * sixpi2_1_3 * (rhomot .* rhomot .* rho_updn);
rho_inv = rhomot .* rhomot .* rhomot;
coeffss = (1.0/4.0) * sixpi2m1_3 * sixpi2m1_3 * (rho_inv .* rho_inv .* rhomot .* rhomot);
ss = (sigma/4.0) .* coeffss;

if iflag == 1 || iflag == 2 || iflag == 4
    divss = 1.0./(1.0 + mu_divkappa * ss);
    dfxdss = mu * (divss .* divss);
elseif iflag == 3
    divss = exp(-mu_divkappa * ss);
    dfxdss = mu * divss;
end

fx = 1.0 + kappa * (1.0 - divss);
dssdn = (-8.0/3.0) * (ss .* rho_inv);
dfxdn = dfxdss .* dssdn;
dssdg = 2.0 * coeffss;
dfxdg = dfxdss .* dssdg;

ex = ex_lsd .* fx;
v1x = ex_lsd .* ((4.0/3.0) * fx + rho_updn .* dfxdn);
v2x = 0.5 * ex_lsd .* rho_updn .* dfxdg;
end

function [ec,v1c,v2c] = pbec(rho,sigma,iflag)
% pbe correlation 
% @param rho  = total electron density
% @param sigma = |\nabla rho|^2
% iflag=1  J.P.Perdew, K.Burke, M.Ernzerhof, PRL 77, 3865 (1996)
% iflag=2  PBEsol: J.P.Perdew et al., PRL 100, 136406 (2008)
% iflag=3  RPBE: B. Hammer, et al., Phys. Rev. B 59, 7413 (1999)
assert(iflag == 1 || iflag == 2 || iflag == 3);

% parameter 
beta_ = [0.066725 0.046 0.066725];
beta = beta_(iflag);
rsfac = 0.6203504908994000; % (0.75/pi)^(1/3)
sq_rsfac = sqrt(rsfac);
sq_rsfac_inv = 1.0/sq_rsfac;
third = 1.0/3.0;
twom1_3 = 2.0^(-third);
ec0_aa = 0.031091; 
ec0_a1 = 0.21370;  
ec0_b1 = 7.5957;
ec0_b2 = 3.5876;   
ec0_b3 = 1.6382;   
ec0_b4 = 0.49294;  
gamma = (1.0 - log(2.0)) /pi^2;
gamma_inv = 1/gamma;
coeff_tt = 1.0/(4.0 * 4.0 / pi * (3.0 * pi^2)^third);

% computation
rho_updn = rho/2.0;
rho_updnm1_3 = rho_updn.^(-third);
rhom1_3 = twom1_3 * rho_updnm1_3;

rhotot_inv = rhom1_3 .* rhom1_3 .* rhom1_3;
rhotmo6 = sqrt(rhom1_3);
rhoto6 = rho .* rhom1_3 .* rhom1_3 .* rhotmo6;

rs = rsfac * rhom1_3;
sqr_rs = sq_rsfac * rhotmo6;
rsm1_2 = sq_rsfac_inv * rhoto6;

%        Formulas A6-A8 of PW92LSD
ec0_q0 = -2.0 * ec0_aa * (1.0 + ec0_a1 * rs);
ec0_q1 = 2.0 * ec0_aa *(ec0_b1 * sqr_rs + ec0_b2 * rs + ec0_b3 * rs .* sqr_rs + ec0_b4 * rs .* rs);
ec0_q1p = ec0_aa * (ec0_b1 * rsm1_2 + 2.0 * ec0_b2 + 3.0 * ec0_b3 * sqr_rs + 4.0 * ec0_b4 * rs);
ec0_den = 1.0./(ec0_q1 .* ec0_q1 + ec0_q1);
ec0_log = -log(ec0_q1 .* ec0_q1 .* ec0_den);
ecrs0 = ec0_q0 .* ec0_log;
decrs0_drs = -2.0 * ec0_aa * ec0_a1 * ec0_log - ec0_q0 .* ec0_q1p .* ec0_den;

%        Add LSD correlation functional to GGA exchange functional
ec = ecrs0;
v1c = ecrs0 - (rs/3.0) .* decrs0_drs;

%        -----------------------------------------------------------------------------
%        Eventually add the GGA correlation part of the PBE functional
%        Note : the computation of the potential in the spin-unpolarized
%        case could be optimized much further. Other optimizations are left to do.

%        From ec to bb
bb = ecrs0 * gamma_inv;
dbb_drs = decrs0_drs * gamma_inv;

%        From bb to cc
exp_pbe = exp(-bb);
cc = 1.0./(exp_pbe - 1.0);
dcc_dbb = cc .* cc .* exp_pbe;
dcc_drs = dcc_dbb .* dbb_drs;

%        From cc to aa
coeff_aa = beta * gamma_inv;
aa = coeff_aa * cc;
daa_drs = coeff_aa * dcc_drs;

%        Introduce tt : do not assume that the spin-dependent gradients are collinear
dtt_dg = 2.0 * rhotot_inv .* rhotot_inv .* rhom1_3 * coeff_tt;
%        Note that tt is (the t variable of PBE divided by phi) squared
tt = 0.5 * sigma .* dtt_dg;

%        Get xx from aa and tt
xx = aa .* tt;
dxx_drs = daa_drs .* tt;
dxx_dtt = aa;

%        From xx to pade
pade_den = 1.0./(1.0 + xx .* (1.0 + xx));
pade = (1.0 + xx) .* pade_den;
dpade_dxx = -xx .* (2.0 + xx) .* (pade_den.^2);
dpade_drs = dpade_dxx .* dxx_drs;
dpade_dtt = dpade_dxx .* dxx_dtt;

%        From pade to qq
qq = tt .* pade;
dqq_drs = tt .* dpade_drs;
dqq_dtt = pade + tt .* dpade_dtt;

%        From qq to rr
arg_rr = 1.0 + beta * gamma_inv * qq;
div_rr = 1.0./arg_rr;
rr = gamma * log(arg_rr);
drr_dqq = beta * div_rr;
drr_drs = drr_dqq .* dqq_drs;
drr_dtt = drr_dqq .* dqq_dtt;

%        The GGA correlation energy is added
ec = ec + rr;

%        From hh to the derivative of the energy wrt the density
drhohh_drho = rr - third * rs .* drr_drs - (7.0/3.0) * tt .* drr_dtt; %- zeta * dhh_dzeta 
v1c = v1c + drhohh_drho;

%        From hh to the derivative of the energy wrt to the gradient of the
%        density, divided by the gradient of the density
%        (The v3.3 definition includes the division by the norm of the gradient)

v2c = rho .* dtt_dg .* drr_dtt;
end
