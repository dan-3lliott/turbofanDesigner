function [Po3, To3, Wdot_comp] = compressor(givens, Po3f, To3f, Pr_comp)
    %=====ASSUMPTIONS=====
    etaPoly = 0.9; %given in prompt
    Cp = givens.R * (givens.gamma_comp/(givens.gamma_comp-1));

    %=====CALCULATIONS=====
    %outlet Po
    Po3 = Po3f * Pr_comp;
    %adiabatic efficiency
    num = (Pr_comp)^((givens.gamma_comp - 1)/givens.gamma_comp) - 1;
    den = (Pr_comp)^((givens.gamma_comp - 1)/(etaPoly*givens.gamma_comp)) - 1;
    etaAdia = num/den;
    %outlet To
    To3s = To3f * Pr_comp^((givens.gamma_comp - 1)/givens.gamma_comp);
    To3 = ((To3s + To3f*(etaAdia - 1))/etaAdia);
    %Wdot
    Wdot_comp = Cp*(To3 - To3f);
end