% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generates a SOFT distribution function from a combination of an
% analytical avalanche distribution and a quadratic radial profile.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETER VALUES
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
OUTPUTFILE_SOFT = 'soft_test_distribution.mat';
OUTPUTFILE_LUKE = 'luke_test_distribution.mat';
EHat = 4;
Zeff = 2;
lnLambda = 17;

rmin = 1.6;
rmax = 2.2;
pmin = 0.1;
pmax = 100;

nr = 40;
np = 100;
nxi = 80;

m = 9.10938356e-31;
c = 299792458;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINE FUNCTIONAL FORMS
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
g = @(p) sqrt(p.^2 + 1);
A = @(p) (EHat + 1) / (Zeff + 1) .* g(p);
g0 = lnLambda * sqrt(Zeff + 5);

fAva = @(p,xi) m*c*A(p) ./ (2*pi*g0.*p.^2) .* exp(-g(p)/g0 - A(p).*(1-abs(xi))) ./ (1 - exp(-2.0*A(p)));
fr = @(r) (1 - (r-rmin)/(rmax-rmin));

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
r = linspace(rmin, rmax, nr);
p = linspace(pmin, pmax, np);
xi = linspace(-1, 1, nxi);
punits = 'normalized';

[P,XI] = meshgrid(p, xi);

f = zeros(np*nxi, nr);
lf = zeros(np, nxi, nr);
for i=1:nr
    tf = fAva(P,XI)';
    f(:,i) = fr(r(i)) * reshape(tf, [np*nxi,1]);
    lf(:,:,i) = fr(r(i)) * tf;
end
fp0 = fAva(p', xi(1));
fxi0 = fAva(p(1), xi');
fr0 = fr(r') * fp0(1);

if sum(abs(fr0 - f(1,1:nr))) > nr*eps
    error('Hash error: r');
end
if sum(abs(fp0 - f(1:np,1))) > np*eps
    error('Hash error: p');
end
if sum(abs(fxi0 - f(1:np:(nxi*np),1))) > np*eps
    error('Hash error: xi');
end

save(OUTPUTFILE_SOFT, '-v7.3', 'r', 'p', 'xi', 'f', 'fr0', 'fp0', 'fxi0', 'punits');

%% Generate LUKE distribution function
betath_ref = rand;
ne_ref = 1;

f = lf;
mhu = xi;
xrhoG = (r-rmin) / (rmax-rmin);
pn = p/betath_ref;

save(OUTPUTFILE_LUKE, '-v7.3', 'betath_ref', 'f', 'mhu', 'xrhoG', 'pn');
