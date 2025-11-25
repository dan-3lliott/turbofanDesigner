function [Po2, To2] = diffuser(givens)
    %=====ASSUMPTIONS=====
    etaAdia = 0.92;

    %=====CALCULATIONS=====
    if (givens.M >= 1)
        rd = 1 - 0.075*(givens.M - 1)^1.35;
    else
        rd = 1;
    end
    To2 = givens.Ta*(1+((givens.gamma_diff-1)/2)*givens.M^2);
    Po2 = rd*givens.Pa*(1+(etaAdia*(givens.gamma_diff-1)/2)*givens.M^2)^(givens.gamma_diff/(givens.gamma_diff-1));
end