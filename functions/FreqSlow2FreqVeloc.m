function [ F_int,V_int,R_Tv ] = FreqSlow2FreqVeloc( f,p,R_fp )
    % Transform radon panel from frequency and slowness to period and velocity,
    % by 2D interpolation to new grid
    %
    
    [F, P] = meshgrid(f,p);
    
    % interpolate freq and velocity
    V = 1./P;
    v = 1./p;
    % Find finest period and velocity to interpolate
    dv   = min(abs(diff(v)));
    f_int = f;
    v_int   = [min(v):dv:max(v)];
    [F_int, V_int] = meshgrid(f_int,v_int);
    R_Tv = interp2(F,V,R_fp,F_int,V_int);
    
    end
    
    