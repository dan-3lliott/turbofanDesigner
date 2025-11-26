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
givens.eta_comb = 0.99; %combustion efficiency in main combustor
givens.eta_ab = 0.96; %combustion efficiency in afterburner

%define other constraints
requirements.To4_max = 1300; %K
requirements.To6_max = 2200; %K
requirements.b_max = 0.12;
requirements.Pr_comp_max = 60;
requirements.Pr_fan_min = 1.1;
requirements.Pr_fan_max = 1.5;
requirements.beta_max = 12; %this is from the slides, project never gave

%=====DEBUGGING - NO OPTIMIZATION=====
engineStruct = analyzeEngine([30 1.2 2.0 0.018 0.010 0.1], givens);
disp('Testcase Engine:');
disp(engineStruct);

%=====OPTIMIZATION=====
%define guesses for design variables
Pr_comp_guess = 15;
Pr_fan_guess = 1.4;
beta_guess = 2;
f_guess = 0.01;
f_ab_guess = 0.01;
b_guess = 0.1;
x0 = [Pr_comp_guess, Pr_fan_guess, beta_guess, f_guess, f_ab_guess, b_guess];

%define upper and lower bounds for design variables
lb = [0 requirements.Pr_fan_min 0 0 0 0];
ub = [requirements.Pr_comp_max requirements.Pr_fan_max requirements.beta_max inf inf requirements.b_max];

%=====CRUISE ENGINE=====

%define cruise flight condition inputs for optimization
givens.Ta = 220; %K
givens.Pa = 29e3; %Pa
givens.M = 0.86;
requirements.ST = 0.86e3; %N*s/kg

%begin optimization
options = optimoptions("fmincon",...
"Algorithm","interior-point",...
"EnableFeasibilityMode",true,...
"SubproblemAlgorithm","cg",...
"MaxFunctionEvaluations",1e6,...
"MaxIterations",1e6);
objectiveFunction = @(x) calculateTSFC(x, givens);
nonlcon = @(x) constraints(x, givens, requirements);
[out, tsfc, ~, ~, ~, ~, ~] = fmincon(objectiveFunction, x0, [], [], [], [], lb, ub, nonlcon, options);

%outputs
engineStruct = analyzeEngine(out, givens);

%calculate max fuel-air ratios
fmax = sym('fmax');
To4_max_adjusted = requirements.To4_max + givens.Cb1*sqrt(engineStruct.b/requirements.b_max);
Cp_comb = givens.R*(givens.gamma_comb/(givens.gamma_comb-1));
To4_expression = ((1-engineStruct.b)*engineStruct.To3 + (fmax*givens.eta_comb*givens.Q_fuel/Cp_comb))/(1+fmax-engineStruct.b);
implicitSol = To4_expression == To4_max_adjusted;
engineStruct.f_max = double(vpasolve(implicitSol, fmax, [0 inf]));

Cp_ab = givens.R*(givens.gamma_ab/(givens.gamma_ab-1));
To6_expression = ((1+engineStruct.f_max)*engineStruct.To5m + (fmax*givens.eta_ab*givens.Q_fuel/Cp_ab))/(1+engineStruct.f_max+fmax);
implicitSol = To6_expression == requirements.To6_max;
engineStruct.f_ab_max = double(vpasolve(implicitSol, fmax, [0 inf]));

%save final optimized result
disp('Cruise-Optimized Engine:');
disp(engineStruct);
writetable(struct2table(engineStruct), 'data outputs/cruise_optimized_engine.xlsx');

%objective function
function tsfc = calculateTSFC(x, givens)
    %perform analysis
    engineStruct = analyzeEngine(x, givens);

    %tsfc objective
    tsfc = engineStruct.TSFC;
end

%nonlinear constraint function
function [c, ceq] = constraints(x, givens, requirements)
    %perform analysis
    engineStruct = analyzeEngine(x, givens);
    %calculate specific thrust and constrain
    ceq(1) = requirements.ST - engineStruct.ST;
    %calculate To4 and constrain
    To4_max_adjusted = requirements.To4_max + givens.Cb1*sqrt(engineStruct.b/requirements.b_max);
    c(1) = engineStruct.To4 - To4_max_adjusted;
    %constrain To6
    c(2) = engineStruct.To6 - requirements.To6_max;
end

%unpacking function
function [Pr_comp, Pr_fan, beta, f, f_ab, b] = unpack(x)
    Pr_comp = x(1);
    Pr_fan = x(2);
    beta = x(3);
    f = x(4);
    f_ab = x(5);
    b = x(6);
end

%engine analysis function
function engineStruct = analyzeEngine(x, givens)
    %=====COMPONENT FRONT-TO-BACK ANALYSIS=====
    
    %unpack design variables
    [engineStruct.Pr_comp, engineStruct.Pr_fan, engineStruct.beta, engineStruct.f, engineStruct.f_ab, engineStruct.b] = unpack(x);
    %diffuser exit properties
    [engineStruct.Po2, engineStruct.To2] = diffuser(givens);
    %fan exit properties
    [engineStruct.Po3f, engineStruct.To3f, engineStruct.Wdot_fan, engineStruct.d] = fan(givens, engineStruct.Po2, engineStruct.To2, engineStruct.Pr_fan, engineStruct.beta);
    %fan nozzle exit properties
    [engineStruct.uef, engineStruct.Tef, engineStruct.Pef] = fan_nozzle(givens, engineStruct.Po3f, engineStruct.To3f);
    %compressor exit properties
    [engineStruct.Po3, engineStruct.To3, engineStruct.Wdot_comp] = compressor(givens, engineStruct.Po3f, engineStruct.To3f, engineStruct.Pr_comp);
    %fuel pump work
    [engineStruct.Wdot_fpump, engineStruct.Pf2] = fpump(givens, engineStruct.Po3, engineStruct.f, engineStruct.f_ab);
    %combustor exit properties
    [engineStruct.Po4, engineStruct.To4] = combustor(givens, engineStruct.Po3, engineStruct.To3, engineStruct.f, engineStruct.b);
    %turbine exit properties - no bleed
    [engineStruct.Po51, engineStruct.To51] = turbine(givens, engineStruct.Wdot_comp + engineStruct.Wdot_fpump, engineStruct.To4, engineStruct.Po4, engineStruct.b, engineStruct.f);
    %mixed turbine exit properties
    [engineStruct.Po5m, engineStruct.To5m] = turbine_mixer(givens, engineStruct.Po51, engineStruct.To51, engineStruct.b, engineStruct.Po3, engineStruct.To3);
    %fan turbine
    [engineStruct.Po52, engineStruct.To52] = turbine(givens, engineStruct.Wdot_fan, engineStruct.To5m, engineStruct.Po5m, 0, engineStruct.f);
    %afterburner exit properties
    [engineStruct.Po6, engineStruct.To6] = afterburner(givens, engineStruct.Po52, engineStruct.To52, engineStruct.f, engineStruct.f_ab);
    %nozzle exit properties and velocity
    [engineStruct.ue, engineStruct.Te, engineStruct.Pe] = nozzle(givens, engineStruct.Po6, engineStruct.To6);
    
    %=====MISCELLANEOUS PARAMETERS=====

    %calculate vehicle velocity
    engineStruct.u = givens.M*sqrt(givens.gamma_diff*givens.R*givens.Ta);
    
    %calculate specific thrust and fuel consumption
    engineStruct.ST = ((1+engineStruct.f+engineStruct.f_ab)*engineStruct.ue - engineStruct.u) + engineStruct.beta*(engineStruct.uef - engineStruct.u) - engineStruct.d; %assuming pe = pa and applying knockdown from bypass fan
    engineStruct.TSFC = (engineStruct.f+engineStruct.f_ab)/engineStruct.ST;

    %calculate efficiencies
    engineStruct.etath = (((1 + engineStruct.f + engineStruct.f_ab) * (((engineStruct.ue)^2) / 2) - (((engineStruct.u)^2) / 2)) + engineStruct.beta*(engineStruct.uef^2 - engineStruct.u^2)/2) / ((engineStruct.f + engineStruct.f_ab) * givens.Q_fuel);
    num = engineStruct.ST*engineStruct.u;
    term1 = (1+engineStruct.f+engineStruct.f_ab)*(engineStruct.ue^2 - engineStruct.u^2)/2;
    term2 = engineStruct.beta*(engineStruct.uef^2 - engineStruct.u^2)/2;
    den = term1+term2;
    engineStruct.etap = (num/den);
    engineStruct.etao = engineStruct.etap * engineStruct.etath;
end
