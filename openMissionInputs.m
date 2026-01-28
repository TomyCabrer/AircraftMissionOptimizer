function [ ...
    Wo, Wp, We, Wf_max, ...
    M_cruise1, M_cruise2, ...
    R_alt, E_loiter, h_loiter, ...
    n_full, n_loiter, ...
    mean_h1_guess, mean_h2_guess, ...
    CL1_guess, CL2_guess, ...
    M_loi_1_guess, M_loi_2_guess ] = unpackMissionOpts(mission, opts)
%UNPACKMISSIONOPTS  Unpack mission + solver option structs into variables.
%
% Usage:
%   [Wo, Wp, We, Wf_max, M1, M2, R_alt, E_loiter, h_loiter, ...
%    n_full, n_loiter, h1g, h2g, CL1g, CL2g, Mloi1g, Mloi2g] = unpackMissionOpts(mission, opts);

% -------------------- mission --------------------
Wo       = mission.Wo;
Wp       = mission.Wp;
We       = mission.We;
Wf_max   = mission.Wf_max;

M_cruise1 = mission.M_cruise1;
M_cruise2 = mission.M_cruise2;

R_alt    = mission.R_alt;
E_loiter = mission.E_loiter;
h_loiter = mission.h_loiter;

% -------------------- opts --------------------
n_full   = opts.n_full;
n_loiter = opts.n_loiter;

mean_h1_guess = opts.mean_h1_guess;
mean_h2_guess = opts.mean_h2_guess;

CL1_guess = opts.CL1_guess;
CL2_guess = opts.CL2_guess;

M_loi_1_guess = opts.M_loi_1_guess;
M_loi_2_guess = opts.M_loi_2_guess;

end