%% CONJUGATE GRADIENT STEP - JBR 7/12/17
%
% Performs single conjugate gradient (CG) step
%
% Claerbout 1992, pg. 142 (CG step)
% Ji 2006, algorithm 5 (CGG with residual and model weights guide)
%
%
% 
% VARIABLES: d  [Y] - data
%            r  [R] - residual
%            dr [G] - conjugate gradient vector
%            ds [S] - descent vector
%            m  [x] - model solution
%            dm [g] - gradient vector
%            Dm [s] - previous descent step
%
% J. Russell
% github.com/jbrussell

function [m,r,Dm,ds] = cgstep(niter,m,dm,Dm,r,dr,ds)

if niter == 0
    Dm = zeros(size(m));
    ds = zeros(size(r));
    if r'*r == 0
        error('r = 0');
    end
    alpha = (dr'*r)/(dr'*dr);
    beta = 0;
else
    dr_dr = dr'*dr;
    ds_ds = ds'*ds;
    dr_ds = dr'*ds;
    determ = dr_dr*ds_ds - dr_ds*dr_ds + (0.00001*(dr_dr*ds_ds)+1e-15);
    dr_r = dr'*r;
    ds_r = ds'*r;
    alpha = (ds_ds*dr_r - dr_ds*ds_r)/determ;
    beta = (-dr_ds*dr_r + dr_dr*ds_r)/determ;
end

% Model step (Dm)
Dm = alpha*dm + beta*Dm;

% Conjugate (ds)
ds = alpha*dr + beta*ds;

% Update Solution
m = m + Dm;

% Update Residual
r = r - ds;




end
