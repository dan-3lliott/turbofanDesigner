clear all
clc
%=====INPUTS=====

%define testcase flight condition inputs
givens.Ta = 220; %K
givens.Pa = 10e3; %Pa
givens.M = 1.5;
requirements.ST = 2.8e3; %N*s/kg

%define other givens
givens.Pf1 = 104e3; %Pa
givens.rho_fuel = 780; %kg/m^3
givens.dP_finj = 550e3; %Pa
givens.MW = 28.8; %g/mol
givens.R = 8314/givens.MW; %J/kgK
givens.gamma_diff = 1.4;
givens.gamma_fan = 1.4;
givens.gamma_comp = 1.38;
givens.gamma_comb = 1.33;
givens.gamma_turb = 1.33;
givens.gamma_fanturb = 1.33;
givens.gamma_nozz = 1.35;
givens.gamma_fannozz = 1.4;
givens.gamma_ab = 1.32;
givens.Q_fuel = 45e6; %J/kg
givens.Pr_b = 0.98; %burner pressure ratio
givens.Cb1 = 700; %K
givens.Cbeta1 = 0.245e3; %N*s/kg
givens.Patm = 101325; %Pa
givens.Pr_ab = 0.97;

%define other constraints
requirements.To4_max = 1300; %K
requirements.To6_max = 2200; %K
requirements.b_max = 0.12;
requirements.Pr_comp_max = 60;
requirements.Pr_fan_min = 1.1;
requirements.Pr_fan_max = 1.5;
requirements.beta_max = 5; %this is from the slides, project never gave

%=====DEBUGGING - NO OPTIMIZATION=====
[Po3, To3, Wdot_comp] = compressor(givens, 39.34e3, 338, 30)