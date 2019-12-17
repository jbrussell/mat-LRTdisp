function [ PER_int,V_int,R_Tv ] = FreqSlow2PeriodVeloc( f,p,R_fp )
% Transform radon panel from frequency and slowness to period and velocity,
% by 2D interpolation to new grid
%

[F, P] = meshgrid(f,p);

% interpolate period and velocity
PER = 1./F;
per = 1./f;
V = 1./P;
v = 1./p;
% Find finest period and velocity to interpolate
dper = min(abs(diff(per)));
dv   = min(abs(diff(v)));
per_int = [min(per):dper:max(per)];
v_int   = [min(v):dv:max(v)];
[PER_int, V_int] = meshgrid(per_int,v_int);
R_Tv = interp2(PER,V,R_fp,PER_int,V_int);

end

