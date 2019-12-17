%% CONJUGATE GRADIENT METHOD - JBR 7/12/17
%
% Solves conjugate gradient using cgstep() by minimizing sum(r^2)
%
% Claerbout 1992, pg. 143 (cgmeth)
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

function [m,r,niter,fiter] = CGsimple(m0,LL,d,maxiter,rthresh)

m = zeros(size(m0));
%m = m0;
r = d;
Dm = [];
ds = [];

niter = 0;
while niter < maxiter && norm(r)/norm(d) >= rthresh
    dm = LL'*r;
    dr = LL*dm;
    [m,r,Dm,ds] = cgstep(niter,m,dm,Dm,r,dr,ds);
    niter = niter + 1;
end

if niter == maxiter
    fiter = 1; %did not reach minimum residual threshold
elseif niter < maxiter
    fiter = 0; %reached minimum residual threshold before maxiter
elseif niter > maxiter
    fiter = -1; %somehow went beyond maxiter...
end


end
