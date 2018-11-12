% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generates a CODE distribution function by first running CODE for
% a simple primary generation scenario. Aside from the main CODE
% distribution file, a separate file containing values of the
% distribution function evaluated in a number of random points is
% created.
%
% The following scripts must be in the MATLAB path:
%
%   runaway/adam/CODE/CODESettings.m
%   runaway/adam/CODE/CODE_timeDependent.m
%   runaway/adam/scripts/helper_scripts/LegendrePolynomials.m
%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc; clear;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETER VALUES
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
OUT_CORR = 'code_corr_values.mat';
OUT_DIST = 'code_test_distribution.mat';

Ny  = 3000;
Nxi = 120;

timeIndex = 5;
Nevals = 300;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% RUN CODE
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
T = 1500;   % eV
n = 5e19;   % m^{-3}
E = 0.55;   % V/m
Z = 1;
yMax = 120;
tMax = 800;
dt = tMax / 99;

o = CODESettings();
o.SetPhysicalParameters(T,n,Z,E,0,0,0,0,0);
o.SetResolutionParameters(Nxi, Ny, yMax, dt, tMax, 'normalized');
o.sourceMode = 0;
o.stepSkip = 10;
o.nStepsToReturn = 5;

s = CODE_timeDependent(o);

xlim([-yMax, yMax]);
ylim([1e-40,1]);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% EVALUATE DISTRIBUTION
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
pMax = s.delta * s.yMax;
f = zeros(Nevals,1);
p = rand(Nevals,1) * pMax;
% xi < 0 are unreliable, so we only consider xi > 0
%xi = rand(Nevals,1)*2 - 1;
xi = rand(Nevals,1);
fLMatrix = reshape(s.f(:,timeIndex), [Ny,Nxi]);

% Below is based on implementation in SYRUP
plusMin = ones(1,Nxi);
plusMin(2:2:end) = -1;
feval = @(p,xi) sum((LegendrePolynomials(Nxi-1, xi)' .*...
    interp1((s.y*s.delta)', fLMatrix, p))*diag(plusMin), 2);

for i=1:Nevals
    f(i) = feval(p(i), xi(i));
end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAVE EVERYTHING
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
save(OUT_CORR, 'p', 'xi', 'f', 'timeIndex');
save(OUT_DIST, '-struct', 's');
