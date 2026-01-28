function [LD_cruise, CL_cruise, CL_star, LD_max] = LDMetircSolver(h, M, cp, CL_guess,engine_on,uc_on,flaps)
%LDMetircSolver  Total drag / L-D metrics at (h, M) using constants in cp.
%
% Inputs
%   h        altitude
%   M        Mach number
%   cp       struct of constants (see fields used below)
%   CL_guess initial CL guess for the internal fixed-point iteration
%   engine_on are both engines on (false = only one engine on)
%   uc_on is the undercarrige deployed (true = undercarrige deployed)
%   flaps are the flaps deployed (true = deployed)
%
% Outputs
%   LD_cruise  cruise L/D (defined here as 0.866 * LD_max)
%   CL_cruise  CL corresponding to LD_cruise
%   CL_star    CL at maximum L/D
%   LD_max     maximum L/D
%
% Dependencies: Re_find.m, Cf_find.m, FF_finder.m, miss_drag.m

% -------------------- Unpack (only what is used) -------------------------
f_fuselage  = cp.f_fuselage;
f_motor     = cp.f_motor;
t_c         = cp.t_c;
x_c_max     = cp.x_c_max;
Lambda_max  = deg2rad(cp.Lambda_max);
A_u_c_front = cp.A_u_c_front;
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

delta_f_deg = cp.delta_f;


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


% Basic input checks
if ~isscalar(h) || ~isscalar(M) || ~isscalar(CL_guess)
    error('h, M, and CL_guess must be scalars.');
end

% -------------------- Compressibility factor keM --------------------------
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
CD1 = sum(Cf .* FF .* Q .* Swet) / S_ref;

delta_f = deg2rad(delta_f_deg);

if uc_on 
    A_u_c_front = A_u_c_front;
else
    A_u_c_front = 0;
    
end


[CDmiss, ~, ~, ~, ~] = miss_drag(A_u_c_front, S_ref, b_f, flaps, delta_f, A_EFF, lengths(1), beta, engine_on, cp);

CD0 = CD1 + CDmiss;

% leakage/protuberance correction (as in your original)
CDlp = (0.035*CD0) / (1 - 0.035);
CD0  = CD0 + CDlp;

% -------------------- Fixed-point iteration for CL_cruise -----------------
CL_cruise = CL_guess;

for it = 1:40
    % finite-wing lift curve slope
    CLalpha = 0.98*(2*pi*AR) ./ (2 + sqrt( (AR.^2 .* (1 + tand(phi50).^2 - M.^2)/0.89^2) + 4 ));

    % Oswald-type efficiency model (your form)
    P = (CLalpha./CL_cruise).*theta.*v + ((CLalpha.*theta).^2 ./ (CL_cruise.^2)).*w + 0.38*CD0;
    e = (keM ./ (Q_e + P*pi*AR)) * ewingl;

    % induced factor
    K = 1/(pi*AR*e);

    % extra drag terms (flap + wave)
    
    if flaps 
        CD_else = (k_f^2) * (delta_CL^2) * cos(Sweep_c/4);
    else
        CD_else = 0;
    end
    if M > Mcrit
        CD_else = CD_else + 20*(M - Mcrit)^4 * (Sc/S_ref);
    end

    CD0_eff = CD0 + CD_else;

    % L/D max and CL at L/D max
    CL_star = sqrt(CD0_eff / K);
    LD_max  = 1 / (2*sqrt(K*CD0_eff));

    % cruise L/D (your definition)
    LD_cruise = 0.866 * LD_max;

    % Solve CL / (CD0_eff + K CL^2) = LD_cruise analytically:
    % (K*LD) CL^2 - CL + (LD*CD0_eff) = 0
    a = K * LD_cruise;
    bq = -1;
    c = LD_cruise * CD0_eff;

    disc = bq^2 - 4*a*c;

    if disc < 0
        % fallback: keep previous CL_cruise (shouldn’t happen for sane inputs)
        break;
    end

    r1 = (1 + sqrt(disc)) / (2*a);
    r2 = (1 - sqrt(disc)) / (2*a);

    % pick the positive root closest to CL_star
    candidates = [r1, r2];
    candidates = candidates(candidates > 0 & isfinite(candidates));
    if isempty(candidates)
        break;
    end
    [~, idx] = min(abs(candidates - CL_star));
    CL_new = candidates(idx);

    % update
    CL_cruise = CL_new;
end
end