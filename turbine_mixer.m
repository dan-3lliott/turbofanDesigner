function [Po5m, To5m] = turbine_mixer(givens, Po51, To51, b, Po3, To3)
    %mixed stag temp using bleed air from compressor exit
    To5m = (1-b)*To51 + b*To3;
    %turbine exhaust stag pressure
    %etaMix = 0.14228;
    %Po5m = (1-b)*Po51 + etaMix*b*Po3;
    Po5m = Po51;
end