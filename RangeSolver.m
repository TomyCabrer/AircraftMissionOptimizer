function [hcruise_start, hcruise_end, hcruise2_start, hcruise2_end, ...
          h_cruise, h_cruise2, ...
          W3cruise_vec, W6cruise_vec, W8loiter_vec, ...
          Wf_Wo, Wf_Wo_real, R_main,W7decloit, W5climb] = ...
          RangeSolver(LD1, LD2, M, M2, CL, CL_2, L_Dmax, ...
                     Wp, Wo, We, Wf_max, ...
                      cp, R_alt, E_loiter, h_loiter, mean_h1, mean_h2)
%RangeSolver  Solves mission range by matching fuel fraction constraints.
%
% Outputs NaN if payload + empty weight exceed take-off weight.

% -------------------- Constants & Environment ----------------------------
g0 = 9.80665;
Re = 6.371e6;
g_of_h = @(h) g0 * (Re / (Re + h))^2;

% Laod aircraft data
SFC_cru  = cp.SFC_cru;
SFC_loit = cp.SFC_loi;
S        = cp.S_ref;

SFC_loiter = SFC_loit * g_of_h(h_loiter);

% -------------------- Mission weight fractions (fixed) -------------------
W1warm  = cp.W1warm;
W2climb = cp.W2climb;
W4dec   = cp.W4dec;
W10land = cp.W10land;

% -------------------- Feasibility checks --------------------------------
Wf_Wo_real = (Wo - Wp - We) / Wo;
Wf_Wo_max  = Wf_max / Wo;

if (Wp + We) > Wo
    [hcruise_start, hcruise_end, hcruise2_start, hcruise2_end, ...
     h_cruise, h_cruise2, ...
     W3cruise_vec, W6cruise_vec, W8loiter_vec, ...
     Wf_Wo, Wf_Wo_real, R_main] = deal(NaN);
    return
end

if Wf_Wo_real > Wf_Wo_max || isnan(Wf_Wo_real)
    Wo = Wp + We + Wf_max;
    Wf_Wo_real = Wf_max / Wo;
end

% -------------------- Initial conditions --------------------------------
WstartCru1 = W1warm * W2climb;

R_main = 1e6;   % initial guess [m]
n = 0;

% -------------------- Nested range search --------------------------------
for i = 1:3
    
    scale = 1e6 / (10^i);
    check = -1;
    R_main = R_main - 1e6 / (10^(i-1));
    
    while check < 0
        R_main = R_main + scale;
        
        % iterate to converge internal mission fractions
        for it = 1:10
            [W3cruise, W6cruise, W8loiter, ...
             hcruise_start, hcruise_end, ...
             hcruise2_start, hcruise2_end, ...
             W7decloit, W5climb] = ...
                WeightFracFinder(LD1, LD2, CL, CL_2, L_Dmax, ...
                            Wo, WstartCru1, ...
                            SFC_cru, SFC_loiter, ...
                            M, M2, S, ...
                            R_main, R_alt, E_loiter, cp, mean_h1, mean_h2);
        end
        
        % total fuel fraction used
        Wf_Wo = 1.01 * (1 - ...
                 W1warm * W2climb * W3cruise * W4dec * ...
                 W5climb * W6cruise * W7decloit * ...
                 W8loiter * W10land);
        
        check = Wf_Wo - Wf_Wo_real;
    end
end

% -------------------- Store outputs --------------------------------------
n = n + 1;

W3cruise_vec(n) = W3cruise;
W6cruise_vec(n) = W6cruise;
W8loiter_vec(n) = W8loiter;

h_cruise  = 0.5 * (hcruise_start  + hcruise_end);
h_cruise2 = 0.5 * (hcruise2_start + hcruise2_end);

end