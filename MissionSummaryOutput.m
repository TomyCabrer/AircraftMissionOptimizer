function MissionSummaryOutput(cp, mission, opts, out)
%MISSIONSUMMARYOUTPUT  Print a readable mission summary to the MATLAB console.
%
% Purpose
%   Post-process solver output (`out`) and print:
%     - solved cruise altitudes and range
%     - key aerodynamic metrics (CL, L/Dmax)
%     - segment weight fractions and fuel closure
%     - convergence/consistency diagnostics
%     - small tables for quick inspection
%
% Inputs
%   cp      constants struct (contains fixed Raymer-style fractions, etc.)
%   mission mission definition (weights, Mach numbers, loiter setup, etc.)
%   opts    solver options / initial guesses (not heavily used here)
%   out     solver output struct returned by your mission solver
%
% Notes
%   - This function assumes `unpackMissionOpts` and `unpackOut` exist.
%   - Units:
%       weights in N, altitude in m, distance in m, time in s
%   - `clc` is called to keep output clean; remove if undesirable for users.

clc;

% -------------------- Unpack mission + opts (only what is used) -----------
[ ...
    Wo, ~, ~, Wf_max, ...
    M_cruise1, M_cruise2, ...
    R_alt, E_loiter, h_loiter, ...
    ~, ~, ...
    ~, ~, ...
    ~, ~, ...
    ~, ~] = unpackMissionOpts(mission, opts);

% -------------------- Unpack solver output --------------------------------
[ ...
    idx, R_main, ...
    hcruise_start, hcruise_end, hcruise2_start, hcruise2_end, ...
    ~, ~, ...
    CL1, CL2, ...
    M_loi_1, M_loi_2, ...
    CL3, CL4, ...
    LDmax_loi, LDmax_loi_2, ...
    W1cruise, W2cruise, Wloiter, ...
    a_loi, ~, ~,...
    ~, ~, ~, ~, ...
    ~, ~, ~, ...
    Wf_Wo, Wf_Wo_real, ...
    h1_diff, h2_diff, ...
    M_loi_diff_1, M_loi_diff_2, ...
    LD1diff, LD2diff, LD3diff, LD4diff ] = unpackOut(out);

% -------------------- Derived quantities ----------------------------------
% Loiter true airspeed from Mach and ISA speed of sound at loiter altitude.
V_loi_start = M_loi_1 * a_loi;
V_loi_end   = M_loi_2 * a_loi;

% Fuel fraction closure (best-case index)
fuel_frac_used  = Wf_Wo(idx);
fuel_frac_avail = Wf_Wo_real(idx);
fuel_margin     = fuel_frac_avail - fuel_frac_used;

% Simple mission weight bookkeeping using fixed fractions in cp + solved fractions.
W0                = Wo;
W_after_warm_climb = Wo * cp.W1warm * cp.W2climb;
W_after_cruise1    = W_after_warm_climb * W1cruise;
W_after_dec        = W_after_cruise1 * cp.W4dec;

% Placeholder: if you model missed approach / climb explicitly, replace here.
W_after_climb2     = W_after_dec;

W_after_cruise2    = W_after_climb2 * W2cruise;
W_after_loiter     = W_after_cruise2 * Wloiter;
W_after_land       = W_after_loiter * cp.W10land;

% Fuel used (simple bookkeeping)
fuel_used_N    = Wo - W_after_land;
fuel_used_kg   = fuel_used_N / 9.80665;
fuel_cap_kg    = Wf_max / 9.80665;
fuel_margin_kg = fuel_cap_kg - fuel_used_kg;

t_loiter_min   = E_loiter / 60;

% -------------------- Console summary -------------------------------------
fprintf("============================================================\n");
fprintf("                 MISSION SOLUTION SUMMARY\n");
fprintf("============================================================\n");

fprintf("\n-- Mission inputs --\n");
fprintf("Cruise 1 Mach: %.3f | Cruise 2 Mach: %.3f\n", M_cruise1, M_cruise2);
fprintf("Alternate cruise distance (R_alt): %.0f km\n", R_alt/1000);
fprintf("Loiter: %.0f min at %.0f m\n", t_loiter_min, h_loiter);

fprintf("\n-- Solved cruise altitudes --\n");
fprintf("Cruise 1: start %.0f m -> end %.0f m (mean %.0f m)\n", ...
    hcruise_start, hcruise_end, 0.5*(hcruise_start+hcruise_end));
fprintf("Cruise 2: start %.0f m -> end %.0f m (mean %.0f m)\n", ...
    hcruise2_start, hcruise2_end, 0.5*(hcruise2_start+hcruise2_end));

fprintf("\n-- Aerodynamics --\n");
fprintf("CL cruise1: %.3f | CL cruise2: %.3f\n", CL1, CL2);
fprintf("Loiter CL: start %.3f | end %.3f\n", CL3, CL4);
fprintf("Loiter L/Dmax: start %.3f | end %.3f\n", LDmax_loi, LDmax_loi_2);

fprintf("\n-- Segment weight fractions (W_after/W_before) --\n");
fprintf("Warmup: %.4f | Climb1: %.4f | Cruise1: %.4f | Decel: %.4f | Cruise2: %.4f | Loiter: %.4f | Land: %.4f\n", ...
    cp.W1warm, cp.W2climb, W1cruise, cp.W4dec, W2cruise, Wloiter, cp.W10land);

fprintf("\n-- Range / fuel closure --\n");
fprintf("Main range solved: %.0f km\n", R_main/1000);
fprintf("Fuel fraction used:  %.5f\n", fuel_frac_used);
fprintf("Fuel fraction avail: %.5f\n", fuel_frac_avail);
fprintf("Fuel fraction margin: %.5f\n", fuel_margin);
fprintf("Fuel used: %.0f kg | Fuel cap: %.0f kg | Margin: %.0f kg\n", ...
    fuel_used_kg, fuel_cap_kg, fuel_margin_kg);

fprintf("\n-- Diagnostics --\n");
fprintf("Cruise altitude convergence (last iter): h1 %.3f %% | h2 %.3f %%\n", h1_diff, h2_diff);
fprintf("Loiter Mach convergence (last iter): start %.3f %% | end %.3f %%\n", M_loi_diff_1, M_loi_diff_2);
fprintf("L/D consistency (vs CD_total_f3):\n");
fprintf("  Cruise1 start %.3f %% | end %.3f %%\n", LD1diff, LD2diff);
fprintf("  Cruise2 start %.3f %% | end %.3f %%\n", LD3diff, LD4diff);

fprintf("\nLoiter speed: start %.1f m/s | end %.1f m/s\n", V_loi_start, V_loi_end);
fprintf("============================================================\n");

% -------------------- Tables (quick inspection) ---------------------------
Seg = ["W0"; "After warm+climb1"; "After cruise1"; "After decel"; ...
       "After cruise2"; "After loiter"; "After land"];
W_N  = [W0; W_after_warm_climb; W_after_cruise1; W_after_dec; ...
        W_after_cruise2; W_after_loiter; W_after_land];
W_kg = W_N / 9.80665;

T_weights = table(Seg, W_N, W_kg, 'VariableNames', {'Segment','Weight_N','Weight_kg'});
disp(T_weights);

AltSeg = ["Cruise1 start"; "Cruise1 end"; "Cruise2 start"; "Cruise2 end"; "Loiter"];
Alt_m  = [hcruise_start; hcruise_end; hcruise2_start; hcruise2_end; h_loiter];
T_alt  = table(AltSeg, Alt_m, 'VariableNames', {'Point','Altitude_m'});
disp(T_alt);

end