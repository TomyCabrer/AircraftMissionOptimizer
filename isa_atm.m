function [a,rho,theta]=isa_atm(h)
T0=288.15; p0=101325; R=287.05287; g=9.80665; gamma=1.4; L=-0.0065;
if h<=11000
    T=T0+L*h; p=p0*(T/T0)^(-g/(L*R));
else
    T=216.65; p=p0*( (T0+L*11000)/T0 )^(-g/(L*R)) * exp(-g*(h-11000)/(R*T));
end
rho=p/(R*T); theta=T/288.15; a=sqrt(gamma*R*T);
end