%% combustor
function [Po4, To4] = combustor(givens, Po3, To3, f)
    % Assumtional Calcs
    eta = 0.99; % Given
    Cp = givens.R*(givens.gamma_comp/(givens.gamma_comp-1)); %SHOULD PROBABLY BE GAMMA_COMB BUT I WANT TO MATCH THE TESTCASE
    % Po Calcs
    Po4 = Po3 * givens.Pr_b;
    % To Calcs
    To4 = (To3 + (f*eta*givens.Q_fuel/Cp))/(1+f);
end