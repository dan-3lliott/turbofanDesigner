function [Po6, To6] = afterburner(givens, Po5m, To5m, f_ab)
    %=====ASSUMPTIONS=====
    eta = 0.96;
    Cp = givens.R*(givens.gamma_ab/(givens.gamma_ab-1));

    %=====CALCULATIONS=====
    %exit stagnation pressure
    Po6 = Po5m * givens.Pr_ab;
    %exit stagnation temperature
    To6 = (Cp*To5m + eta * f_ab * givens.Q_fuel) / ((1 + f_ab)*Cp);
end