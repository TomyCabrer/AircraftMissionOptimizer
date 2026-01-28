% List of component names (fuselage, wing, engines, vertical tail, horizontal tail)
clear
cp = make_const_params('Const_Param.mat');

f_fuselage =cp.f_fuselage  ;
f_motor    =cp.f_motor    ;
t_c        =cp.t_c        ;
x_c_max    =cp.x_c_max    ;
Lambda_max =cp.Lambda_max ;
b          =cp.b          ;
b_f        =cp.b_f        ;
A_u_c_front=cp.A_u_c_front;
S_ref      =cp.S_ref      ;



A_EFF      =cp.A_EFF      ;
beta       =cp.beta       ;
Mcrit      =cp.Mcrit      ;
Sc         =cp.Sc         ;
AR         =cp.AR         ;
k_f        =cp.k_f        ;

Sweep_c    =cp.Sweep_c    ;
smooth_surf=cp.smooth_surf;
t_c_vec=cp.t_c_vec;
x_c_max_vec=cp.x_c_max_vec;
perc_lam = cp.perc_lam;
engine_on=false; 
lengths=cp.lengths;
Q=cp.Q;
what_component=cp.what_component;
components_name=cp.components_name;
Swet=cp.Swet;

delta_f    =cp.delta_f_lan    ;
delta_CL = cp.delta_CL;


CL=2.9;
h=9000;
h=0;
M=0.78;
% Landing
%M= 0.1483;
%take off
%M=0.1693;
M=0.1;
Lambda_max=Lambda_max*pi/180;
for i = 1:length(what_component)
    [Re(i,1), ~, ~] = Re_find(h, M, lengths(i), smooth_surf);
    [Cf(i,1)] = Cf_find(Re(i),M,perc_lam(i));
    FF(i,1) = FF_finder(what_component(i), f_fuselage, f_motor, M, t_c_vec(i), x_c_max_vec(i), Lambda_max);
end

ae=cp.ae;
be=cp.be;
ce=cp.ce;
Mcomp=cp.Mcomp;
e_theo=cp.e_theo;
keM = ae .* ((M./Mcomp - 1) .^ be) + ce;
keM(M < Mcomp) = 1;              % below compressibility onset
keM = max(0, min(1, keM));
e=e_theo*keM;

ewingl =cp.ewingl;
Q_e= cp.Q_e;
theta=cp.theta;
v=cp.v;
w=cp.w;
phi50=cp.phi50;
%theta=-2;



% The drag coefficient is normalized by the total wetted surface area
CD1 = sum(Cf .* FF .* Q .* Swet) / S_ref; % Sum of drag coefficients from individual components
% Convert the fuselage upsweep angle (beta) from degrees to radians
beta = beta * pi / 180;
% Calculate the miscellaneous drag contributions (underbody drag, flap drag, etc.)
[CDmiss, C_D_u_c, C_D_HLD, C_D_WE, C_D_upsweep] = miss_drag(A_u_c_front, S_ref, b_f, b, delta_f, A_EFF, lengths(1), beta, engine_on);

% Total drag coefficient (CD0) is the sum of the component-based drag and miscellaneous drag
CD0 = CD1 + CDmiss;
CDlp=(0.035*CD0)/(1-0.035);
CD0=CD0 +CDlp;

CD0=0.0906;

CLalpha = 0.98*(2*pi*AR) ./ (2 + sqrt((AR.^2 .* (1 + (tand(phi50)).^2 - M.^2)/0.89^2+4)));
P = (CLalpha./CL).*theta.*v + ((CLalpha.*theta).^2 ./ CL.^2).*w + 0.38*CD0;
e = keM ./ (Q_e + P*pi*AR)*ewingl;% -cor;

K = 1/(pi * AR *e);
%K_h = 1/(pi * AR_h *e_h);

CD = CD0 + K * CL^2; %+ eta_h * (S_h / S_ref) * K_h * CL_h^2;


% Calculate the change in induced drag (Delta C_Di)
% Formula: Delta C_Di = k_f^2 * (Delta C_L)^2 * cos(Lambda_c / 4)
Sweep_c= Sweep_c*pi/180;
if delta_f ~= 0 % Check if we have falps 
    delta_CD_i = k_f^2 * delta_CL^2 * cos(Sweep_c);
else
    delta_CD_i = 0;
end

CD = CD + delta_CD_i;

if M>Mcrit
 CD = CD + 20*(M - Mcrit)^4*Sc/S_ref;
end

