function [CD0, CD, LD] = DragFinder(h, M, cp, CL, engine_on, uc_on, flaps)
%DragFinder  Total drag coefficient and L/D at (h, M, CL).
%
% Inputs
%   h         altitude [m]
%   M         Mach number [-]
%   cp        constants struct
%   CL        lift coefficient [-]
%   engine_on true: both engines on; false: one engine on
%   uc_on     true: undercarriage deployed; false: retracted
%   flaps     true: flaps deployed; false: retracted
%
% Outputs
%   CD0       zero-lift drag coefficient [-]
%   CD        total drag coefficient [-]
%   LD        lift-to-drag ratio [-]
%
% Dependencies: Re_find.m, Cf_find.m, FF_finder.m, miss_drag.m

% -------------------- Basic checks ---------------------------------------
if ~(isscalar(h) && isscalar(M) && isscalar(CL))
    error('h, M, and CL must be scalars.');
end

% -------------------- Unpack constants -----------------------------------
f_fuselage  = cp.f_fuselage;
f_motor     = cp.f_motor;
t_c         = cp.t_c;
x_c_max     = cp.x_c_max;
Lambda_max  = deg2rad(cp.Lambda_max);

A_u_c_front_base = cp.A_u_c_front;
b_f         = cp.b_f;
S_ref       = cp.S_ref;
A_EFF       = cp.A_EFF;
beta        = deg2rad(cp.beta);

Mcrit       = cp.Mcrit;
Sc          = cp.Sc;
AR          = cp.AR;

k_f         = cp.k_f;
delta_CL    = cp.delta_CL;
Sweep_c     = deg2rad(cp.Sweep_c);
delta_f     = deg2rad(cp.delta_f);   % flap deflection angle used by miss_drag

smooth_surf = cp.smooth_surf;
perc_lam    = cp.perc_lam;
lengths     = cp.lengths;
Q           = cp.Q;
what_comp   = cp.what_component;
Swet        = cp.Swet;

% compressibility efficiency model
ae    = cp.ae;  be = cp.be;  ce = cp.ce;
Mcomp = cp.Mcomp;

% wing efficiency model params
ewingl = cp.ewingl;
Q_e    = cp.Q_e;
theta  = cp.theta;
v      = cp.v;
w      = cp.w;
phi50  = cp.phi50;

% -------------------- Compressibility factor -----------------------------
keM = ae .* ((M./Mcomp - 1) .^ be) + ce;
keM(M < Mcomp) = 1;
keM = min(1, max(0, keM));  % clamp to [0,1]

% -------------------- Per-component Re, Cf, FF ----------------------------
n  = numel(what_comp);
Re = zeros(n,1);
Cf = zeros(n,1);
FF = zeros(n,1);

for i = 1:n
    Re(i) = Re_find(h, M, lengths(i), smooth_surf);
    Cf(i) = Cf_find(Re(i), M, perc_lam(i));
    FF(i) = FF_finder(what_comp(i), f_fuselage, f_motor, M, t_c, x_c_max, Lambda_max);
end

% -------------------- Parasite (zero-lift) drag ---------------------------
CD_parasite = sum(Cf .* FF .* Q .* Swet) / S_ref;

A_u_c_front = A_u_c_front_base;
if ~uc_on
    A_u_c_front = 0;
end

[CDmiss, ~, ~, ~, ~] = miss_drag( ...
    A_u_c_front, S_ref, b_f, flaps, delta_f, A_EFF, lengths(1), beta, engine_on, cp);

CD0 = CD_parasite + CDmiss;

% leakage/protuberance correction (kept from your model)
CD0 = CD0 + (0.035*CD0) / (1 - 0.035);

% -------------------- Induced drag model ---------------------------------
CLalpha = 0.98*(2*pi*AR) ./ (2 + sqrt((AR.^2 .* (1 + tand(phi50).^2 - M.^2)/0.89^2) + 4));

P = (CLalpha./CL).*theta.*v + ((CLalpha.*theta).^2 ./ (CL.^2)).*w + 0.38*CD0;
e = (keM ./ (Q_e + P*pi*AR)) * ewingl;

K = 1 / (pi * AR * e);

CD = CD0 + K * CL.^2;

% -------------------- Incremental drag terms -----------------------------
if flaps
    CD = CD + (k_f^2) * (delta_CL^2) * cos(Sweep_c);
end

if M > Mcrit
    CD = CD + 20*(M - Mcrit)^4 * (Sc/S_ref);
end

LD = CL / CD;
end