function [Po6, To6] = afterburner(givens, Po5m, To5m, f, f_ab)
    %=====ASSUMPTIONS=====
    eta = givens.eta_ab;
    Cp = givens.R*(givens.gamma_ab/(givens.gamma_ab-1));

    %=====CALCULATIONS=====
    %exit stagnation pressure
    Po6 = Po5m * givens.Pr_ab;
    %exit stagnation temperature
    To6 = ((1+f)*To5m + (f_ab*eta*givens.Q_fuel/Cp))/(1+f+f_ab);
end