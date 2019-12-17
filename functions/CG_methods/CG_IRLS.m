%% CONJUGATE GRADIENT METHOD - JBR 7/12/17
%
% Solves conjugate gradient using cgstep() by minimizing sum(r^2)
%
% Claerbout 1992, pg. 143 (cgmeth)
%
% Ji 2006, algorithm 2 [iteratively reweighted least squares (IRLS)]
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

function [m,r,niter,fiter] = CG_IRLS(m0,LL,d,maxiter,rthresh)

m = zeros(size(m0));
%m = m0;
r = d;
Dm = [];
ds = [];

epsilon = max(abs(d))/100; % to avoid divide by zero
% epsilon = prctile(d,2); % to avoid divide by zero

niter = 0;
while niter < maxiter && norm(r)/norm(d) >= rthresh
    
    % Weighting Matrices
    if niter == 0
        Wm = eye(length(m));
        Wr = eye(length(r));
    else
%         Wm = diag(abs(m).^(1.5)); % L-1 norm
%         Wm = diag(abs(m)); % L0 norm
        Wm = diag(abs(m).^(1/2)); % L1 norm
%         Wm = diag(ones(size(m))); % L2 norm
        Wr = diag(abs(r).^(-1/2)); % L1 norm
%         Wr = diag(abs(r).^(-1)); % L0 norm
        Wr(diag(abs(r))<=epsilon & diag(abs(r))~=0) = epsilon; % Fix division by zero in Wr
    end
    r = Wr*(LL*Wm*m-d);
    
    niter_inner = 0;
    while niter_inner < maxiter && norm(r)/norm(d) >= rthresh  
        dm = Wm'*LL'*Wr'*r;
        dr = Wr*LL*Wm*dm;
%         dr = LL*dm;
        [m,r,Dm,ds] = cgstep(niter,m,dm,Dm,r,dr,ds);
        niter_inner = niter_inner + 1;
    end
    m = Wm*m;
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
