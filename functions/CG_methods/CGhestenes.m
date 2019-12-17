%% CONJUGATE GUIDED GRADIENT (HESTENES METHOD) - JBR 10/29/15
%
% Performs conjugate guided gradient (CGG) method using more numerically stable algorithm
%
% Claerbout 1992, pg. 151 (CG Hestenes)
% Ji 2006, algorithm 5 (CGG with residual and model weights guide)
%
% INPUT: m0 - starting model (only the size is used, not the value)
% 		 LL - L matrix
%		  d - data matrix
%	maxiter - maximum iterations to perform
%   rthresh - minimum residual required to break out of function
%
% OUTPUT: m - final model vector
%		  r - final residual vector
%     niter - number of iterations performed
%     fiter - iteration flag : 1 => did not reach minimum residual threshold
%     						   0 => reached rthresh before maxiter
%							  -1 => PROBELM! somehow went beyond maxiter...
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

function [m,r,niter,fiter] = CGhestenes(m0,LL,d,maxiter,rthresh)

m = zeros(size(m0));
%m = m0;
r = d;
dm = LL'*d;
Dm = dm;
gamma_neg = dm'*dm;

% epsilon = max(abs(d))/100; % to avoid divide by zero
epsilon = prctile(d,2); % to avoid divide by zero

niter = 0;
while niter < maxiter && norm(r)/norm(d) >= rthresh
    niter = niter + 1;
    ds = LL*Dm;
    alpha = gamma_neg/(ds'*ds);
    m = m + alpha*Dm;
    r = r - alpha*ds;
    Wm = diag(abs(m).^(1/2)); % L1 norm
%     Wm = diag(ones(size(m))); % L2 norm
    Wr = diag(abs(r).^(-1/2)); % L1 norm
    % Fix division by zero in Wr
    Wr(diag(abs(r))<=epsilon & diag(abs(r))~=0) = epsilon;
    %
    dm  = Wm'*LL'*Wr'*r; % weighted
    %dm = LL'*r;
    gamma = dm'*dm;
    beta = gamma/gamma_neg;
    gamma_neg = gamma;
    Dm = dm + beta*Dm;
end

if niter == maxiter
    fiter = 1; %did not reach minimum residual threshold
elseif niter < maxiter
    fiter = 0; %reached minimum residual threshold before maxiter
elseif niter > maxiter
    fiter = -1; %somehow went beyond maxiter...
end


end
