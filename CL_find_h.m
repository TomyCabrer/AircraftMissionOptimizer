function f = CL_find_h(h2, CL, W2, M, S)
%CL_FIND_H  Residual of the lift balance used to solve cruise altitude.
%
% Solves for the altitude at which the lift equals the aircraft weight:
%   0.5 * rho(h) * V(h)^2 * S * CL = W
% with V = M * a(h).
%
% Used with fzero to find the altitude where steady, level flight is possible.

    % Atmospheric properties at altitude h2
    [a2, rho2, ~] = isa_atm(h2);

    % True airspeed from Mach number
    V2 = M * a2;

    % Lift-to-weight ratio (should be 1 at equilibrium)
    lift_ratio = (S * CL * rho2 * V2^2) / (2 * W2);

    % Residual for root finding
    f = lift_ratio - 1;
end