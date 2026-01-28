function [Re,Re1,Re2] = Re_find(h,M,l,K)
%   h  - Altitude (m)
%   M  - Mach number
%   l  - Characteristic length (m)
%   K  - Characteristic dimension (m), used for calculating Re2
% Output:
%   Re  - The maximum Reynolds number (either Re1 or Re2)
%   Re1 - Reynolds number based on Mach and air density
%   Re2 - Reynolds number based on a critical value and characteristic length

[a,rho,theta]=isa_atm(h);
mu0 = 1.7894e-5;  % Dynamic viscosity at sea level (kg/m·s)
T0 = 288.15;  % Reference temperature (K)
C = 198.72;  % Sutherland's constant (K)

% Calculate dynamic viscosity using Sutherland's law
mu = mu0 * (theta)^1.5 * ((T0 + C) / (theta*T0 + C));
V=M.*a;
Re1=rho.*V.*l./mu;
Re2=38.21*(l./K).^(1.053);

if Re1<Re2
    Re=Re1;
else
    Re=Re2;
end

end