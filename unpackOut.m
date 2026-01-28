function [ ...
    idx, R_main, ...
    hcruise_start, hcruise_end, hcruise2_start, hcruise2_end, ...
    mean_h1, mean_h2, ...
    CL1, CL2, ...
    M_loi_1, M_loi_2, ...
    CL3, CL4, ...
    LDmax_loi, LDmax_loi_2, ...
    W1cruise, W2cruise, Wloiter, ...
    a_loi, W7decloit, W5climb,...
    h1s, h1e, h2s, h2e, ...
    W1_vec, W2_vec, Wloi_vec, ...
    Wf_Wo, Wf_Wo_real, ...
    h1_diff, h2_diff, ...
    M_loi_diff_1, M_loi_diff_2, ...
    LD1diff, LD2diff, LD3diff, LD4diff ] = unpackOut(out)
%UNPACKOUT  Convenience function to unpack the solver output struct.
%
% This function extracts all relevant quantities stored in the `out` struct
% (returned by `runMissionSolver`) into individual variables. It is intended
% for users who prefer working with flat variables instead of nested structs.
%
% Typical usage:
%   out = runMissionSolver(cp, mission, cfg, opts);
%   [idx, R_main, hcruise_start, hcruise_end, ...] = unpackOut(out);
%
% This is especially useful for:
%   - Post-processing
%   - Plotting
%   - Parametric studies
%   - Exporting results
%
% The `out` struct is assumed to have the following fields:
%   out.solution     → main scalar results for the selected best case
%   out.vectors      → full vectors across all candidate cases
%   out.diagnostics  → convergence and consistency metrics
%
% The outputs are grouped as follows:
%
% 1) Best-case scalars (from out.solution)
%    - idx              : index of the selected best case
%    - R_main           : solved main cruise range [m]
%    - hcruise_start    : start altitude of cruise 1 [m]
%    - hcruise_end      : end altitude of cruise 1 [m]
%    - hcruise2_start   : start altitude of cruise 2 [m]
%    - hcruise2_end     : end altitude of cruise 2 [m]
%    - mean_h1, mean_h2 : mean altitudes of cruise segments [m]
%    - CL1, CL2         : lift coefficients at cruise 1 and cruise 2 [-]
%    - M_loi_1, M_loi_2 : loiter Mach numbers (start and end) [-]
%    - CL3, CL4         : loiter lift coefficients (start and end) [-]
%    - LDmax_loi        : maximum L/D at loiter (start) [-]
%    - LDmax_loi_2      : maximum L/D at loiter (end) [-]
%    - W1cruise         : weight fraction for cruise 1 [-]
%    - W2cruise         : weight fraction for cruise 2 [-]
%    - Wloiter          : weight fraction for loiter [-]
%    - a_loi            : speed of sound at loiter altitude [m/s]
%    - W7decloit        : weight fraction 2nd decent
%    - W5climb          : weight fraction 2nd climb
%
% 2) Full vectors (from out.vectors)
%    - h1s, h1e         : start/end altitudes for cruise 1 over all cases [m]
%    - h2s, h2e         : start/end altitudes for cruise 2 over all cases [m]
%    - W1_vec           : cruise 1 weight fraction vector [-]
%    - W2_vec           : cruise 2 weight fraction vector [-]
%    - Wloi_vec         : loiter weight fraction vector [-]
%    - Wf_Wo            : fuel fraction required for each case [-]
%    - Wf_Wo_real       : available fuel fraction for each case [-]
%
% 3) Diagnostics (from out.diagnostics)
%    - h1_diff          : convergence error for cruise 1 altitude [%]
%    - h2_diff          : convergence error for cruise 2 altitude [%]
%    - M_loi_diff_1     : convergence error for loiter Mach (start) [%]
%    - M_loi_diff_2     : convergence error for loiter Mach (end) [%]
%    - LD1diff–LD4diff  : L/D consistency checks against optimal values [%]
%
% This function does not modify the input struct and performs no validation;
% it assumes `out` is a valid output of `runMissionSolver`.

% -------------------- solution --------------------
idx            = out.solution.idx;
R_main         = out.solution.R_main;

hcruise_start  = out.solution.hcruise_start;
hcruise_end    = out.solution.hcruise_end;
hcruise2_start = out.solution.hcruise2_start;
hcruise2_end   = out.solution.hcruise2_end;

mean_h1        = out.solution.mean_h1;
mean_h2        = out.solution.mean_h2;

CL1            = out.solution.CL1;
CL2            = out.solution.CL2;

M_loi_1        = out.solution.M_loi_1;
M_loi_2        = out.solution.M_loi_2;

CL3            = out.solution.CL3;
CL4            = out.solution.CL4;

LDmax_loi      = out.solution.LDmax_loi;
LDmax_loi_2    = out.solution.LDmax_loi_2;

W1cruise       = out.solution.W1cruise;
W2cruise       = out.solution.W2cruise;
Wloiter        = out.solution.Wloiter;

a_loi          = out.solution.a_loi;

W7decloit      = out.solution.W7decloit;
W5climb        = out.solution.W5climb;

% -------------------- vectors --------------------
h1s            = out.vectors.h1s;
h1e            = out.vectors.h1e;
h2s            = out.vectors.h2s;
h2e            = out.vectors.h2e;

W1_vec         = out.vectors.W1_vec;
W2_vec         = out.vectors.W2_vec;
Wloi_vec       = out.vectors.Wloi_vec;

Wf_Wo          = out.vectors.Wf_Wo;
Wf_Wo_real     = out.vectors.Wf_Wo_real;

% -------------------- diagnostics --------------------
h1_diff        = out.diagnostics.h1_diff_pct;
h2_diff        = out.diagnostics.h2_diff_pct;

M_loi_diff_1   = out.diagnostics.M_loi_diff_start_pct;
M_loi_diff_2   = out.diagnostics.M_loi_diff_end_pct;

LD_diffs       = out.diagnostics.LD_diffs_pct;
LD1diff        = LD_diffs(1);
LD2diff        = LD_diffs(2);
LD3diff        = LD_diffs(3);
LD4diff        = LD_diffs(4);

end