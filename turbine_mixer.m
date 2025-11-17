function [Po5m, To5m] = turbine_mixer(Po51, To51, b, To3)
    %mixed stag temp using bleed air from compressor exit
    To5m = (1-b)*To51 + b*To3;
    %turbine exhaust stag pressure is the same
    Po5m = Po51;
end