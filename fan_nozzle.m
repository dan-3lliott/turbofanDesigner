function [uef, Tef, Pef] = fan_nozzle(givens, Po3f, To3f)
% Nozzle (fan). Computes exhaust velocity, exit temperature, and exit pressure.
% Assumptions: adiabatic and isentropic. eta_noz = 0.95.

%% Input/Givens
eta_noz = 0.97;

Pa = givens.Pa;
gamma = givens.gamma_fannozz;
R = givens.R;
Cp = gamma*R/(gamma - 1);
Pef = Pa;

%% Calculations

% Isentropic exit temp
Te_is = To3f * (Pef/Po3f)^((gamma - 1)/gamma);

% Actual Exit velocity
uef = sqrt(eta_noz) * sqrt(2*Cp*(To3f - Te_is));

% Actual exit  temp
Tef = To3f - (uef^2)/(2*Cp);
end
