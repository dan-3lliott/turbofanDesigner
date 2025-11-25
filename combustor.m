%% combustor
function [Po4, To4] = combustor(givens, Po3, To3, f, b)
    % Assumtional Calcs
    eta = givens.eta_comb; % Given
    Cp = givens.R*(givens.gamma_comb/(givens.gamma_comb-1));
    % Po Calcs
    Po4 = Po3 * givens.Pr_b;
    % To Calcs
    To4 = ((1-b)*To3 + (f*eta*givens.Q_fuel/Cp))/(1+f-b);
end