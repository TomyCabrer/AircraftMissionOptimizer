function out = MissionProfileSolver(cp, mission, cfg, opts)
%MISSIONPROFILESOLVER  Core mission-profile solver.
%
% This function couples:
%   - aerodynamic performance (via DragFinder / LDMetircSolver),
%   - weight-fraction bookkeeping,
%   - cruise altitude convergence,
%   - loiter Mach convergence,
%   - and fuel-range closure
% into a single mission-consistent solution.
%
% Inputs
%   cp      : aircraft constant-parameter struct (geometry, aerodynamics, SFC, etc.)
%   mission : mission definition struct (weights, Mach numbers, loiter, distances)
%   cfg     : configuration flags per segment (engines, gear, flaps)
%   opts    : numerical solver options and initial guesses
%
% Output
%   out : struct with three sub-structs:
%       out.solution    → best-case scalar results
%       out.vectors     → vectors of all candidate solutions
%       out.diagnostics → convergence and consistency metrics
%
% The solver proceeds in four main stages:
%   1) Seed loiter L/Dmax
%   2) Iterate cruise altitudes + fuel-range closure
%   3) Iterate loiter Mach number (start and end)
%   4) Select the best feasible fuel solution and build output

% -------------------------------------------------------------------------
% Unpack mission and solver options into flat variables
% -------------------------------------------------------------------------
[ ...
    Wo, Wp, We, Wf_max, ...
    M_cruise1, M_cruise2, ...
    R_alt, E_loiter, h_loiter, ...
    n_full, n_loiter, ...
    mean_h1, mean_h2, ...
    CL1, CL2, ...
    M_loi_1, M_loi_2] = unpackMissionOpts(mission, opts);

% -------------------------------------------------------------------------
% 1) Initial loiter L/Dmax seed (used in the weight-fraction model)
% -------------------------------------------------------------------------
[~, ~, ~, LDmax_loi] = LDMetircSolver(h_loiter, 0.30, cp, CL1, ...
    cfg.loiter.engine_on, cfg.loiter.uc_on, cfg.loiter.flaps);

% -------------------------------------------------------------------------
% 2) Main coupled solve: cruise altitudes + range closure
%    Outer loop enforces consistency between:
%      - aerodynamic performance
%      - weight fractions
%      - cruise altitudes
% -------------------------------------------------------------------------
for k = 1:n_full
    mean_h1_prev = mean_h1;
    mean_h2_prev = mean_h2;

    % Cruise 1 aerodynamics
    [LD1, CL1, ~, ~] = LDMetircSolver(mean_h1, M_cruise1, cp, CL1, ...
        cfg.cruise1.engine_on, cfg.cruise1.uc_on, cfg.cruise1.flaps);

    % Cruise 2 aerodynamics
    [LD2, CL2, ~, ~] = LDMetircSolver(mean_h2, M_cruise2, cp, CL2, ...
        cfg.cruise2.engine_on, cfg.cruise2.uc_on, cfg.cruise2.flaps);

    % Weight fractions + altitude consistency + range closure
    [h1s, h1e, h2s, h2e, mean_h1, mean_h2, ...
     W1_vec, W2_vec, Wloi_vec, Wf_Wo, Wf_Wo_real, R_main, W7decloit, W5climb] = ...
        RangeSolver(LD1, LD2, M_cruise1, M_cruise2, CL1, CL2, LDmax_loi, ...
                   Wp, Wo, We, Wf_max, cp, ...
                   R_alt, E_loiter, h_loiter, mean_h1, mean_h2);
end

% Convergence indicators for cruise altitudes
h1_diff = 100 * (mean_h1_prev - mean_h1) / mean_h1;
h2_diff = 100 * (mean_h2_prev - mean_h2) / mean_h2;

% -------------------------------------------------------------------------
% 3) Loiter Mach number solve (start of loiter)
%    Uses simple lift equilibrium: L = W
% -------------------------------------------------------------------------
[a_loi, ~, ~] = isa_atm(h_loiter);
S = cp.S_ref;  % Reference wing area

for k = 1:n_loiter
    M_prev = M_loi_1;

    [~, ~, CL3, LDmax_loi] = LDMetircSolver(h_loiter, M_loi_1, cp, CL1, ...
        cfg.loiter.engine_on, cfg.loiter.uc_on, cfg.loiter.flaps);

    % Original formulation retained: rho = 1.225, a = 340
    M_loi_1 = sqrt(2 * Wo * W2_vec(end) / (CL3 * 1.225 * S)) / 340;
end
M_loi_diff_1 = 100 * (M_prev - M_loi_1) / M_loi_1;

% -------------------------------------------------------------------------
% 4) Loiter Mach number solve (end of loiter)
% -------------------------------------------------------------------------
for k = 1:n_loiter
    M_prev = M_loi_2;

    [~, ~, CL4, LDmax_loi_2] = LDMetircSolver(h_loiter, M_loi_2, cp, CL1, ...
        cfg.loiter.engine_on, cfg.loiter.uc_on, cfg.loiter.flaps);

    M_loi_2 = sqrt(2 * Wo * Wloi_vec(end) / (CL4 * 1.225 * S)) / 340;
end
M_loi_diff_2 = 100 * (M_prev - M_loi_2) / M_loi_2;

% -------------------------------------------------------------------------
% 5) Select best feasible solution
%    Criterion: minimum mismatch between required and available fuel fraction
% -------------------------------------------------------------------------
[~, idx] = min(abs(Wf_Wo - Wf_Wo_real));

hcruise_start  = h1s(idx);  hcruise_end  = h1e(idx);
hcruise2_start = h2s(idx);  hcruise2_end = h2e(idx);

W1cruise = W1_vec(idx);
W2cruise = W2_vec(idx);
Wloiter  = Wloi_vec(idx);

% -------------------------------------------------------------------------
% 6) Optional L/D consistency checks
%    Compare DragFinder (direct) vs LDMetircSolver (optimised) results
% -------------------------------------------------------------------------
[LDopt1, ~, ~, ~] = LDMetircSolver(hcruise_start, M_cruise1, cp, CL1, ...
    cfg.cruise1.engine_on, cfg.cruise1.uc_on, cfg.cruise1.flaps);
[~, ~, LD1a] = DragFinder(hcruise_start, M_cruise1, cp, CL1, ...
    cfg.cruise1.engine_on, cfg.cruise1.uc_on, cfg.cruise1.flaps);

[LDopt2, ~, ~, ~] = LDMetircSolver(hcruise_end, M_cruise1, cp, CL1, ...
    cfg.cruise1.engine_on, cfg.cruise1.uc_on, cfg.cruise1.flaps);
[~, ~, LD1b] = DragFinder(hcruise_end, M_cruise1, cp, CL1, ...
    cfg.cruise1.engine_on, cfg.cruise1.uc_on, cfg.cruise1.flaps);

[LDopt3, ~, ~, ~] = LDMetircSolver(hcruise2_start, M_cruise2, cp, CL2, ...
    cfg.cruise2.engine_on, cfg.cruise2.uc_on, cfg.cruise2.flaps);
[~, ~, LD2a] = DragFinder(hcruise2_start, M_cruise2, cp, CL2, ...
    cfg.cruise2.engine_on, cfg.cruise2.uc_on, cfg.cruise2.flaps);

[LDopt4, ~, ~, ~] = LDMetircSolver(hcruise2_end, M_cruise2, cp, CL2, ...
    cfg.cruise2.engine_on, cfg.cruise2.uc_on, cfg.cruise2.flaps);
[~, ~, LD2b] = DragFinder(hcruise2_end, M_cruise2, cp, CL2, ...
    cfg.cruise2.engine_on, cfg.cruise2.uc_on, cfg.cruise2.flaps);

LD1diff = 100 * (LD1a - LDopt1) / LDopt1;
LD2diff = 100 * (LD1b - LDopt2) / LDopt2;
LD3diff = 100 * (LD2a - LDopt3) / LDopt3;
LD4diff = 100 * (LD2b - LDopt4) / LDopt4;

% -------------------------------------------------------------------------
% 7) Build structured output
% -------------------------------------------------------------------------
out = struct();

% Best-case scalars
out.solution = struct( ...
    "idx", idx, ...
    "R_main", R_main, ...
    "hcruise_start", hcruise_start, "hcruise_end", hcruise_end, ...
    "hcruise2_start", hcruise2_start, "hcruise2_end", hcruise2_end, ...
    "mean_h1", mean_h1, "mean_h2", mean_h2, ...
    "CL1", CL1, "CL2", CL2, ...
    "M_loi_1", M_loi_1, "M_loi_2", M_loi_2, ...
    "CL3", CL3, "CL4", CL4, ...
    "LDmax_loi", LDmax_loi, "LDmax_loi_2", LDmax_loi_2, ...
    "W1cruise", W1cruise, "W2cruise", W2cruise, "Wloiter", Wloiter, ...
    "a_loi", a_loi, "W7decloit", W7decloit, "W5climb", W5climb);

% Full solution vectors
out.vectors = struct( ...
    "h1s", h1s, "h1e", h1e, "h2s", h2s, "h2e", h2e, ...
    "W1_vec", W1_vec, "W2_vec", W2_vec, "Wloi_vec", Wloi_vec, ...
    "Wf_Wo", Wf_Wo, "Wf_Wo_real", Wf_Wo_real);

% Convergence + consistency diagnostics
out.diagnostics = struct( ...
    "h1_diff_pct", h1_diff, ...
    "h2_diff_pct", h2_diff, ...
    "M_loi_diff_start_pct", M_loi_diff_1, ...
    "M_loi_diff_end_pct", M_loi_diff_2, ...
    "LD_diffs_pct", [LD1diff, LD2diff, LD3diff, LD4diff]);

end