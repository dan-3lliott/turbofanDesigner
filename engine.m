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
requirements.beta_max = 12; %this is from the slides, project never gave

%=====DEBUGGING - NO OPTIMIZATION=====
engineStruct = analyzeEngine([30 1.2 2.0 0.018 0.010 0.1], givens);
disp(engineStruct);

%define actual flight condition inputs for optimization
givens.Ta = 220; %K
givens.Pa = 29e3; %Pa
givens.M = 0.86;
requirements.ST = 0.86e3; %N*s/kg

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
disp(engineStruct);
writetable(struct2table(engineStruct), 'engine_out.xlsx');

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
    engineStruct.Wdot_fpump = fpump(givens, engineStruct.Po3, engineStruct.f, engineStruct.f_ab);
    %combustor exit properties
    [engineStruct.Po4, engineStruct.To4] = combustor(givens, engineStruct.Po3, engineStruct.To3, engineStruct.f);
    %turbine exit properties - no bleed
    [engineStruct.Po51, engineStruct.To51] = turbine(givens, engineStruct.Wdot_fan, engineStruct.Wdot_comp, engineStruct.Wdot_fpump, engineStruct.To4, engineStruct.Po4, engineStruct.b);
    %mixed turbine exit properties
    [engineStruct.Po5m, engineStruct.To5m] = turbine_mixer(engineStruct.Po51, engineStruct.To51, engineStruct.b, engineStruct.To3);
    %afterburner exit properties
    [engineStruct.Po6, engineStruct.To6] = afterburner(givens, engineStruct.Po5m, engineStruct.To5m, engineStruct.f_ab);
    %nozzle exit properties and velocity
    [engineStruct.ue, engineStruct.Te, engineStruct.Pe] = nozzle(givens, engineStruct.Po6, engineStruct.To6);
    
    %=====MISCELLANEOUS PARAMETERS=====

    %calculate vehicle velocity
    engineStruct.u = givens.M*sqrt(givens.gamma_diff*givens.R*givens.Ta);
    
    %calculate specific thrust and fuel consumption
    engineStruct.ST = ((1+engineStruct.f)*engineStruct.ue - engineStruct.u) + engineStruct.beta*(engineStruct.uef - engineStruct.u) - engineStruct.d; %assuming pe = pa and applying knockdown from bypass fan
    engineStruct.TSFC = (engineStruct.f)/engineStruct.ST;

    %calculate efficiencies
    engineStruct.etath = (((1 + engineStruct.f) * (((engineStruct.ue)^2) / 2) - (((engineStruct.u)^2) / 2)) + engineStruct.beta*(engineStruct.uef^2 - engineStruct.u^2)/2) / ((engineStruct.f) * givens.Q_fuel);
    engineStruct.etap = 2 * (((1+engineStruct.f)*(engineStruct.ue/engineStruct.u) - 1)/((1+engineStruct.f)*((engineStruct.ue/engineStruct.u)^2) - 1)); 
    engineStruct.etao = engineStruct.etap * engineStruct.etath;
end
