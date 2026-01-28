%% Mission profile solver – main execution script
% This script defines the mission, numerical options, and aircraft
% configuration, then runs the mission-profile solver and post-processing.
% The only blocks a user normally needs to modify are:
%   - Mission inputs
%   - Solver options
%   - Configuration flags

clear;

% -------------------- Load aircraft constants -----------------------------
% Loads all fixed geometric, aerodynamic, and propulsion parameters.
cp = make_const_params("Const_Param.mat");


%% ===================== USER-DEFINED MISSION INPUTS =====================
% These parameters describe the aircraft mass breakdown and the mission
% performance requirements. Units must be respected.

mission = struct();
opts    = struct();

% ---- Aircraft weights ---------------------------------------------------
mission.Wo     = NaN;        % [N] Maximum take-off weight (MTOW)
mission.Wp     = NaN;        % [N] Payload weight
mission.We     = NaN;        % [N] Operating empty weight
mission.Wf_max = NaN;    % [N] Maximum usable fuel weight

% ---- Cruise conditions --------------------------------------------------
mission.M_cruise1 = NaN;         % [-] Mach number for main cruise
mission.M_cruise2 = NaN;         % [-] Mach number for alternate cruise

% ---- Mission geometry ---------------------------------------------------
mission.R_alt    = NaN;        % [m] Distance to alternate airport, cruise 2 distance
mission.E_loiter = NaN;       % [s] Loiter endurance (45 minutes)
mission.h_loiter = NaN;          % [m] Loiter altitude


%% ===================== SOLVER NUMERICAL OPTIONS =====================
% These control convergence behaviour and provide initial guesses.
% They do NOT constrain the final solution, only how fast and robustly it
% converges.

opts.n_full   = 5;                % Outer iterations (cruise altitude + fuel closure)
opts.n_loiter = 7;                % Iterations for loiter Mach convergence

% ---- Initial guesses ----------------------------------------------------
opts.mean_h1_guess = 9000;        % [m] Initial guess: mean altitude of cruise 1
opts.mean_h2_guess = 7000;        % [m] Initial guess: mean altitude of cruise 2
opts.CL1_guess     = 0.365;       % [-] Initial CL for cruise 1
opts.CL2_guess     = 0.31;        % [-] Initial CL for cruise 2
opts.M_loi_1_guess = 0.30;        % [-] Initial Mach guess at start of loiter
opts.M_loi_2_guess = 0.30;        % [-] Initial Mach guess at end of loiter


%% ===================== AIRCRAFT CONFIGURATION PER SEGMENT =====================
% These flags define the aerodynamic configuration used in each flight phase.
% They directly affect drag and fuel burn.

cfg.cruise1 = struct( ...
    "engine_on", true, ...   % true → all engines operating
    "uc_on",     false, ...  % true → landing gear deployed
    "flaps",     false );    % true → flaps deployed

cfg.cruise2 = struct( ...
    "engine_on", true, ...
    "uc_on",     false, ...
    "flaps",     false );

cfg.loiter = struct( ...
    "engine_on", true, ...
    "uc_on",     false, ...
    "flaps",     false );


%% ===================== SOLVER EXECUTION =====================
% Runs the coupled mission solver:
%   - cruise altitude determination
%   - fuel fraction closure
%   - loiter Mach convergence
%   - consistency checks

out = MissionProfileSolver(cp, mission, cfg, opts);


%% ===================== POST-PROCESSING =====================
% Console report with all main performance metrics
MissionSummaryOutput(cp, mission, opts, out)

% Mission profile plot:
%   - altitude vs mission node
%   - weight fraction vs mission node
MissionPLot(cp, mission, opts, out)
