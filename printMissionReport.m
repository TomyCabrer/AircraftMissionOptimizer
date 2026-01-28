
% ========================= 3) OUTPUT/REPORT FUNCTION =========================
function printMissionReport(results)
%PRINTMISSIONREPORT  Pretty console output + key diagnostics.

sol = results.solution;
diag = results.diagnostics;
m = results.inputs.mission;

fprintf("============================================================\n");
fprintf("                 MISSION SOLUTION SUMMARY\n");
fprintf("============================================================\n");

fprintf("\n-- Inputs --\n");
fprintf("Wo=%.0f N | Wp=%.0f N | We=%.0f N | Wf_max=%.0f N\n", m.Wo, m.Wp, m.We, m.Wf_max);
fprintf("M_cruise1=%.3f | M_cruise2=%.3f\n", m.M_cruise1, m.M_cruise2);
fprintf("R_main=%.0f km | R_alt=%.0f km | Loiter=%.0f min @ %.0f m\n", ...
    m.R_main_m/1000, m.R_alt_m/1000, m.t_loiter_s/60, m.h_loiter_m);

fprintf("\n-- Solved altitudes --\n");
fprintf("Cruise1: start %.0f m -> end %.0f m\n", sol.hcruise_start, sol.hcruise_end);
fprintf("Cruise2: start %.0f m -> end %.0f m\n", sol.hcruise2_start, sol.hcruise2_end);

fprintf("\n-- Solved fractions --\n");
fprintf("W1cruise=%.4f | W2cruise=%.4f | Wloiter=%.4f\n", sol.W1cruise, sol.W2cruise, sol.Wloiter);
fprintf("Range solved R_main = %.0f km\n", sol.R_main/1000);

fprintf("\n-- Loiter --\n");
fprintf("M_loiter start=%.3f | end=%.3f\n", sol.M_loi_1, sol.M_loi_2);
fprintf("CL_loiter start=%.3f | end=%.3f | LDmax start=%.3f | end=%.3f\n", ...
    sol.CL3, sol.CL4, sol.LDmax_loi, sol.LDmax_loi_2);

fprintf("\n-- Diagnostics --\n");
fprintf("Cruise mean height change: h1 %.3f %% | h2 %.3f %%\n", diag.h1_diff_pct, diag.h2_diff_pct);
fprintf("Loiter Mach change: start %.3f %% | end %.3f %%\n", diag.M_loi_diff_start_pct, diag.M_loi_diff_end_pct);
fprintf("L/D diffs (%%): [%.3f %.3f %.3f %.3f]\n", diag.LD_diffs_pct);

fprintf("============================================================\n");
end


% ========================= helpers (small) =========================
function mission = applyMissionDefaults(mission)
if ~isfield(mission,"R_alt_m"), mission.R_alt_m = 370e3; end
if ~isfield(mission,"t_loiter_s"), mission.t_loiter_s = 45*60; end
if ~isfield(mission,"h_loiter_m"), mission.h_loiter_m = 5000*0.3048; end
end

function opts = applyOptsDefaults(opts)
if ~isfield(opts,"n_full"), opts.n_full = 5; end
if ~isfield(opts,"n_loiter"), opts.n_loiter = 7; end
if ~isfield(opts,"M_loi_1_init"), opts.M_loi_1_init = 0.30; end
if ~isfield(opts,"M_loi_2_init"), opts.M_loi_2_init = 0.30; end
end

function validateMissionInputs(mission, cp, cfg, opts) %#ok<INUSD>
mustBePositive(mission.Wo); mustBePositive(mission.Wp); mustBePositive(mission.We);
if mission.Wo <= (mission.Wp + mission.We)
    error("Wo must exceed Wp + We.");
end
mustBeInRange(mission.M_cruise1, 0.1, 0.95);
mustBeInRange(mission.M_cruise2, 0.1, 0.95);
mustBePositive(mission.R_alt_m);
mustBePositive(mission.t_loiter_s);
mustBeNonnegative(mission.h_loiter_m);
if ~isfield(cp,"S_ref") || ~isscalar(cp.S_ref) || cp.S_ref<=0
    error("cp.S_ref must be a positive scalar.");
end
end

function mustBeInRange(x, a, b)
if ~(isscalar(x) && x>=a && x<=b)
    error("Value must be in range [%.3g, %.3g].", a, b);
end
end