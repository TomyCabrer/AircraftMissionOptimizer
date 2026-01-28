function [W3cruise, W6cruise, W8loiter, ...
          hcruise_start, hcruise_end, hcruise2_start, hcruise2_end, ...
          W7decloit, W5climb] = ...
          WeightFracFinder(LD1, LD2, CL, CL_2, L_Dmax, ...
                      Wo, WstartCru1, ...
                      SFC_cru, SFC_loiter, ...
                      M, M2, S, ...
                      R_main, R_alt, E_loiter, cp, ...
                      h_cruise, h_cruise2)
%WEIGHTFRACFINDER  Computes segment weight fractions and cruise altitudes.
%
% Notes:
% - Requires: CL_find_h.m, isa_atm.m
% - Uses g_of_h from cp if provided; otherwise assumes constant g0.

% -------------------- Defaults / constants -------------------------------
W3cruise = cp.W3cruise;   % initial guesses / seeds
W6cruise = cp.W6cruise;
W8loiter = cp.W8loiter;

% gravity model (prefer cp.g_of_h if you store it; else constant g0)
if isfield(cp, "g_of_h") && isa(cp.g_of_h, "function_handle")
    g_of_h = cp.g_of_h;
else
    g0 = 9.80665;
    g_of_h = @(~) g0;
end

% -------------------- Cruise 1 start altitude ----------------------------
fun = @(h) CL_find_h(h, CL, Wo*WstartCru1, M, S);
hcruise_start = fzero(fun, h_cruise);

% -------------------- Iterate: update end weights/altitudes + fractions ---
for it = 1:10

    % cruise 1 end altitude (depends on W3cruise)
    WendCru1 = WstartCru1 * W3cruise;
    fun = @(h) CL_find_h(h, CL, Wo*WendCru1, M, S);
    hcruise_end = fzero(fun, h_cruise);

    % update mean altitude for next fzero guesses
    h_cruise = 0.5 * (hcruise_start + hcruise_end);

    % climb/descend factors between cruises (your original model)
    W5climb   = 1 - (1 - 0.985) * (h_cruise2 / h_cruise);
    W7decloit = 1 - (1 - 0.995) * (h_cruise2 / h_cruise);

    % cruise 2 start altitude
    WstartCru2 = WendCru1 * W7decloit * W5climb;
    fun = @(h) CL_find_h(h, CL_2, Wo*WstartCru2, M2, S);
    hcruise2_start = fzero(fun, h_cruise2);

    % cruise 2 end altitude
    WendCru2 = WstartCru2 * W6cruise;
    fun = @(h) CL_find_h(h, CL_2, Wo*WendCru2, M2, S);
    hcruise2_end = fzero(fun, h_cruise2);

    % update mean altitude for next fzero guesses
    h_cruise2 = 0.5 * (hcruise2_start + hcruise2_end);

    % -------------------- Speeds, times, SFC -----------------------------
    [a1, ~, ~] = isa_atm(h_cruise);
    [a2, ~, ~] = isa_atm(h_cruise2);

    V1 = a1 * M;
    V2 = a2 * M2;

    t1 = R_main / V1;
    t2 = R_alt  / V2;

    SFC1 = SFC_cru * g_of_h(h_cruise);
    SFC2 = SFC_cru * g_of_h(h_cruise2);

    % -------------------- Weight fractions -------------------------------
    W3cruise = exp(-t1 * SFC1 / LD1);
    W6cruise = exp(-t2 * SFC2 / LD2);
    W8loiter = exp(-E_loiter * SFC_loiter / L_Dmax);
end

end