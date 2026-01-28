function [CDmiss, C_D_u_c, C_D_HLD, C_D_WE, C_D_upsweep] = miss_drag(A_u_c_front, S_ref, b_f, flaps, delta_f, A_EFF, D, beta, engine_on, cp)
    % A_u_c_front: Front area of the undercarriage (m²)
    % S_ref: Reference wing area (m²)
    % b_f: Flap span (m)
    % b: Wing span (m)
    % delta_f: Flap deflection angle (degrees)
    % A_EFF: Effective engine area (m²)
    % D: Fuselage diameter (m)
    % beta: Fuselage upsweep angle (dimensionless)

    % Undercarriage drag coefficient (C_D_u_c)
    if A_u_c_front~=0

        n_wheel    = cp.n_wheel;
        area_wheel = cp.area_wheel;
        CD_wheel   = cp.CD_wheel;
        wheel = CD_wheel*(area_wheel / S_ref)*n_wheel;

        n_nwheel    = cp.n_nwheel;
        area_nwheel = cp.area_nwheel;
        CD_nwheel   = cp.CD_nwheel;
        nwheel = CD_nwheel*(area_nwheel / S_ref)*n_nwheel;

        n_strut     =cp.n_strut;
        area_strut  =cp.area_strut;
        CD_strut     =cp.CD_strut;
        strut = CD_strut*(area_strut/S_ref)*n_strut;

        n_nstrut     =cp.n_nstrut;
        area_nstrut  =cp.area_nstrut;
        CD_nstrut     =cp.CD_nstrut;
        nstrut = CD_nstrut*(area_nstrut/S_ref)*n_nstrut;

        
        C_D_u_c = nstrut+strut+nwheel+wheel;
    else
        C_D_u_c = 0;
    end
    
    if flaps
        % Flap drag coefficient (C_D_HLD)
        C_D_HLD = 0.0023 * (b_f) * delta_f;
    else  
        C_D_HLD = 0;
    end
    
    % Windmilling engine drag coefficient (C_D_WE)
    if engine_on
        C_D_WE = 0;
    else
        C_D_WE = 0.3 * (A_EFF / S_ref);
    end
    
    % Fuselage upsweep drag coefficient (C_D_upsweep)
    C_D_upsweep = 3.83 * ((pi * D^2) / (4 * S_ref)) * beta^2.5;
    
    % Total miscellaneous drag coefficient (CDmiss)
    CDmiss = C_D_upsweep + C_D_WE + C_D_HLD + C_D_u_c;
end