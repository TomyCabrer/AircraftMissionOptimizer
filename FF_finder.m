function [FF] = FF_finder(component,f_fuselage,f_motor,M,t_c,x_c_max,Lambda_max)
% what: 1= wing, 2 = fuselage, 3 = nacelle

if component == 1
   
    % Calculate the form factor for the wing 
    FF = (1 + (0.6 .* t_c ./ (x_c_max) + 100 .* (t_c).^4)) ./ (1.34.* M.^0.18 .* cos(Lambda_max).^0.28);

elseif component == 2
    % Calculate the form factor for the fuselage 
    FF = 1 + (60 ./ f_fuselage.^3) + (f_fuselage ./ 400);

elseif component == 3
    % Calculate the form factor for the nacelle 
    FF = 1 + (0.35 ./ f_motor);

else
    error('Invalid component. Choose "Wing", "Fuselage", or "Nacelle".');
end
end