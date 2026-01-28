function [Cf] = Cf_find(Re,M,perc_lam)
  % Cf_find calculates the skin friction coefficient (Cf) based on the Reynolds number (Re),
    % Mach number (M), and the critical Reynolds number (Re_crit).
    %
    % Inputs:
    %   Re      - Reynolds number
    %   M       - Mach number
    %   Re_crit - Critical Reynolds number for transition from laminar to turbulent flow
    %
    % Output:
    %   Cf      - Skin friction coefficient (Cf)

    % If Reynolds number is greater than the critical value, use the turbulent flow formula

    Cf_tur = 0.455 ./ (log10(Re).^2.58 * (1 + 0.144 * M.^2).^0.65);

    Cf_lam=1.326./sqrt(Re);

    Cf=perc_lam * Cf_lam +(1-perc_lam)*Cf_tur;


end