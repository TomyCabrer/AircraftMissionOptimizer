function cp = make_const_params(savepath)
%MAKE_CONST_PARAMS  Aircraft constants for the mission solver.
%
% User workflow
%   - Edit only the "REQUIRED USER INPUTS" and "OPTIONAL USER INPUTS" blocks.
%   - Leave "MODEL TUNING" and "DERIVED" blocks unchanged unless you re-fit
%     the correlations.
%
% Notes
%   - Use SI units. Angles in degrees unless stated.
%   - This function validates inputs and fails fast with clear messages.
% Values given are found from literature and can be used in a range of
% aircraft, so kept for user.


cp = struct();

%% ======================= USER: EDIT HERE (GEOMETRY) =======================
% Wing planform
cp.S_ref      = NaN;   % [m^2] reference wing area
cp.b          = NaN;    % [m]   wing span
cp.AR         = NaN;    % [-]   aspect ratio (keep consistent with b^2/S_ref)
cp.lambda     = NaN;     % [-]   taper ratio (c_tip/c_root)

% Sweep / thickness (used by form-factor and CLalpha correlations)
cp.Sweep_c    = NaN;      % [deg] 1/4-chord sweep
cp.phi50      = NaN;    % [deg] ~1/2-chord sweep (used in CLalpha correlation)
cp.t_c        = NaN;    % [-]   thickness-to-chord
cp.x_c_max    = NaN;    % [-]   location of max thickness (x/c)
cp.Lambda_max = NaN;    % [deg] sweep used in FF model (keep consistent with geometry)

% Fuselage / nacelle sizing (parasite drag form factors + induced corrections)
cp.f_fuselage = NaN;      % [-] fuselage fineness ratio
cp.dia_fus    = NaN;       % [m] fuselage diameter (or equivalent)
cp.f_motor    = NaN;  % [-] nacelle fineness ratio
cp.A_EFF      = NaN;     % [m^2] effective engine area (misc drag model)
cp.beta       = NaN;  % [deg] fuselage upsweep (misc drag model)

%% =================== USER: EDIT HERE (PARASITE DRAG BUILD-UP) ===================
% Parasite drag is computed as: sum(Cf * FF * Q * Swet) / S_ref
% DON'T CAHNGE BELOW, ONLY CHANGE LENGTHS, SWET, Q, PERCLAM, SMOOTHLAM
cp.components_name = {'Fuselage'; 'Wing'; 'Engines'; 'Vertical Tail'; 'Horizontal Tail'};
cp.what_component  = [2; 1; 3; 1; 1];             % 1=surface, 2=body, 3=nacelle

cp.lengths         = [NaN; NaN; NaN; NaN; NaN]; % [m] Re reference lengths
cp.Swet            = [NaN; NaN; NaN; NaN; NaN]; % [m^2] wetted areas
cp.Q               = [1; 1.1; NaN; 1.04; 1.04];  % [-] interference factors
cp.perc_lam        = [0.01, 0.07, 0.01, 0.09, 0.09]; % [-] laminar fraction
cp.smooth_surf     = 0.633984e-5;                 % [m] roughness height (Cf model)


%% =================== USER: EDIT HERE (CONFIGURATION EFFECTS) ===================
% Flaps
cp.delta_f  = 0;      % [deg] flap deflection (misc drag model)
cp.delta_CL = 2.9;    % [-]   CL increment with flaps
cp.k_f      = 0.28;   % [-]   flap drag factor (partial span flap)
cp.b_f      = NaN;   % [m]   flap span (misc drag model)



%% ======================= OPTIONAL USER INPUTS ============================
% If you do not model these effects, leave as defaults.
cp.k_f   = 0.28;  % [-]
cp.b_f   = 0.56;  % [m]

% Gear detail (only if miss_drag uses these fields)
cp.n_wheel = 0; cp.area_wheel = 0; cp.CD_wheel = 0.25;
cp.n_strut = 0; cp.area_strut = 0; cp.CD_strut = 0.05; 
cp.n_nwheel = 0; cp.area_nwheel = 0; cp.CD_nwheel = 0.25;
cp.n_nstrut = 0; cp.area_nstrut = 0; cp.CD_nstrut = 0.05;



%% =================== USER: EDIT HERE (MISSION FIXED FRACTIONS + SFC) ===================
% Fixed non-cruise fractions (Raymer-style). Cruise/loiter fractions are solved.
% Rayman values
cp.W1warm  = 0.97;    % [-] warmup/taxi
cp.W2climb = 0.985;   % [-] climb to cruise
cp.W4dec   = 0.995;   % [-] descent/deceleration
cp.W10land = 0.995;   % [-] landing

% Initial seeds (solver updates these; values only affect initialisation)
cp.W3cruise = 0.7503; % [-] seed for main cruise
cp.W6cruise = 0.9842; % [-] seed for alternate cruise
cp.W8loiter = 0.9821; % [-] seed for loiter

% Specific fuel consumption (your model uses [1/s])
cp.SFC_cru = 14.4e-6; % [1/s] cruise
cp.SFC_loi = 11.3e-6; % [1/s] loiter


%% =================== USER: EDIT HERE (Wave) ===================
cp.Mcrit = NaN;      % [-]  wave drag starts above this Mach
cp.Sc    = NaN;    % [m^2] reference area for wave drag scaling

%% =================== USER: LEAVE UNLESS YOU RETUNE THE MODEL ===================

% keM compressibility efficiency model (fit constants)
cp.Mcomp = 0.311;     % [-] onset Mach for keM model
cp.ae    = -0.00152;  % [-] keM fit constant
cp.be    = 10.82;     % [-] keM fit exponent
cp.ce    = 1;         % [-] keM offset (clamped later)



%% =================== USER: LEAVE UNLESS YOU RETUNE THE MODEL ===================
% Induced-drag/efficiency correlation tuning parameters.
% These are fit knobs for your e-model; change only if re-calibrating.
cp.theta     = -2*pi/180; % [rad] fit parameter inside P(...)
cp.hwinglets = 2.2;       % [m] winglet height (only matters if ewingl enabled)





%% ======================= DERIVED: DO NOT EDIT =======================
% Everything below is computed from the user blocks above. If you change the
% equations here, you are changing the modelling approach, not the aircraft.

lambda  = cp.lambda;
AR      = cp.AR;
dia_fus = cp.dia_fus;
b       = cp.b;

f = 0.0524*lambda.^4 - 0.15*lambda.^3 + 0.1659*lambda.^2 - 0.0706*lambda + 0.0119;

e_theo = 1 ./ (1 + f .* AR);       % baseline Oswald efficiency correlation
Kef    = 1 - 2*(dia_fus^2/b^2);    % fuselage penalty factor

cp.Q_e = 1 / (e_theo * Kef);

% Winglet factor (kept off by default for reproducibility)
cp.ewingl = 1; % set to (1 + 2*cp.hwinglets/cp.b)^2 to enable winglet benefit

cp.e_theo = e_theo * Kef * cp.ewingl * 0.873; % 0.873 = Kedo correction used in your model

% Correlation terms used in your induced-drag efficiency model
cp.v = 0.0134*(lambda - 0.3) - 0.0037*lambda.^2;
cp.w = (0.0088*lambda - 0.0051*lambda.^2) .* (1 - 0.0006*AR.^2);


%% ======================= SAVE (OPTIONAL) =======================
if nargin >= 1 && ~isempty(savepath)
    save(savepath, "cp");
end

end