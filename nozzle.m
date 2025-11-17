function [ue, Te, Pe] = nozzle(givens, Po6, To6)
% Nozzle (core). Computes exhaust velocity, exit temperature, and exit pressure.
% Assumptions: adiabatic and isentropic. eta_noz = 0.95.

%% Input/Givens
eta_noz = 0.95;

Pa = givens.Pa;
gamma = givens.gamma_nozz;
R = givens.R;
Cp = gamma*R/(gamma - 1);
Pe = Pa;

%% Calculations

% Isentropic exit temp
Te_is = To6 * (Pe/Po6)^((gamma - 1)/gamma);

% Actual Exit velocity
ue = sqrt(eta_noz) * sqrt(2*Cp*(To6 - Te_is));

% Actual exit  temp
Te = To6 - (ue^2)/(2*Cp);
end
