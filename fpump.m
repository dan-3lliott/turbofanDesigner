function Wdot_fpump = fpump(givens, Po3, f, f_ab)
    %=====ASSUMPTIONS=====
    etaAdia = 0.35;
    %=====CALCULATIONS=====
    %total pressure ratio
    Pf2 = Po3 + givens.dP_finj;
    dP_fpump = (Pf2 - givens.Pf1);
    Wdot_fpump = ((f+f_ab)*dP_fpump)/(etaAdia*givens.rho_fuel);
end