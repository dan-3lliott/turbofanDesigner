function[Po51, To51] = turbine(givens, Wdot_fan, Wdot_comp, Wdot_fpump, To4, po4, b)
%Note:
%Tmax = 1300K
%b_max = 0.12
%C_bl = 700K


%work reqirement from upstream components 
%Change if wc is not defined to include wp (pump work) 
w_t = Wdot_fan + Wdot_comp + Wdot_fpump; %wdot comp is in J/kg, not a rate
%Accounting for bleed rate from compressor to turbine for cooling
%divinding by 1-b corrects the specific work term

%Givens (Look into structure format from main script for simplification)
gamma_t = givens.gamma_turb; %spec. heat ratio for turbine - from table
Mt = givens.MW; %kg/kmol -  mol.wt turbine gas
R = givens.R; %J/(kg*K) gas constant
eta_pt = 0.92; %polytropic efficiency



%turbine gas property 
cp_t = gamma_t/(gamma_t - 1) * R; %J/(kg*K)

% turbine work = cp_t * (To4 - To5_1) - adiabatic so Qdot is zero
%no change in height and assume V2 =~ V1 so work is change in enthalpy
To51 = To4 - w_t / ((1-b)*cp_t);


%Exit Pressure Calculations - po5
Tri = To51 / To4; % temperature ratio To5 / To4
Pr_t = Tri^(gamma_t / (eta_pt * (gamma_t-1)));  %po5 / po4
Po51 = po4 * Pr_t; % kPa turbine exit total pressure state 5.1
end
