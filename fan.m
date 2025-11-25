function [Po3f, To3f, Wdot_fan, d] = fan(givens, Po2, To2, Pr_fan, beta)
    %=====ASSUMPTIONS=====
    etaPoly = 0.9; %given in prompt
    Cp = givens.R * (givens.gamma_fan/(givens.gamma_fan-1));

    %=====CALCULATIONS=====
    %outlet Po
    Po3f = Po2 * Pr_fan;
    %adiabatic efficiency
    num = (Pr_fan)^((givens.gamma_fan - 1)/givens.gamma_fan) - 1;
    den = (Pr_fan)^((givens.gamma_fan - 1)/(etaPoly*givens.gamma_fan)) - 1;
    etaAdia = num/den;
    %outlet To
    To3fs = To2 * Pr_fan^((givens.gamma_fan - 1)/givens.gamma_fan);
    To3f = ((To3fs + To2*(etaAdia - 1))/etaAdia);
    %Wdot
    Wdot_fan = (1+beta)*Cp*(To3f - To2);
    %drag delta d
    d = (givens.Cbeta1*givens.M^2)*(givens.Pa/givens.Patm)*beta^1.5;
end