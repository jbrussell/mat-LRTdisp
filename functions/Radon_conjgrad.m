function [ Rfft,f ] = Radon_conjgrad(p,t,M,delta,maxiter,rthresh,method)
% Solve radon transform using the conjugate gradient method. 
%
% Options for 4 different conjugate gradient methods:
%
% (1) 'CGsimple': Simple conjugate gradient [*fastest method]
% (2) 'CGG_weight': Conjugate guided gradient (CGG) with residual and model weights guide
% (3) 'CGhestenes': Similar to CGG_weight but uses less precise alternative to cgstep() at each iteration. (* CGG_weight preferred)
% (4) 'CG_IRLS': Iteratively reweighted least squares (IRLS) [*slowest method]
%
% Algorithms from Ji 2006 & Claerbout 1992
% 
% 10/14/19
% J. Russell
% github.com/jbrussell


% Define some array/matrices lengths.
it=length(t);
iF=pow2(nextpow2(it)+1); % Double length
iDelta=length(delta);
ip=length(p);
Mfft=fft(M,iF,2);
dF=1/(t(1)-t(2));
f = ((  [1:floor((iF+1)/2)]  -1)/iF)*dF*-1;

% Define blocks
blk_w = ip; %block width within L
blk_h = iDelta; % block height within L
delta_block = repmat(delta,ip,1)'; %repmat(delta,1,blk_w);
p_block = repmat(p,iDelta,1); %repmat(p',blk_h,1);

m0 = zeros(length(p),1);
Rfft = zeros(ip,iF);
rfft = zeros(size(Mfft));
for j = 1:length(f)
    exp_arg = -1i*2*pi*f(j).*delta_block.*p_block;
    L = exp(exp_arg);
    d = Mfft(:,j);
    switch lower(method) %make it case insensitive, all lower case
        case 'cgsimple'
            [m,r,niter,fiter] = CGsimple(m0,L,d,maxiter,rthresh);
        case 'cgg_weight'
            [m,r,niter,fiter] = CGG_weight(m0,L,d,maxiter,rthresh);
        case 'cghestenes'
            [m,r,niter,fiter] = CGhestenes(m0,L,d,maxiter,rthresh);
        case 'cg_irls'
            [m,r,niter,fiter] = CG_IRLS(m0,L,d,maxiter,rthresh);
        otherwise
            error('Incorrect method name. Choose CGsimple, CGG_weight, CGhestenes, CG_IRLS');
    end
    Rfft(:,j) = m;
    rfft(:,j) = r;
%     Rfft_est(:,j) = L*m; % d_est
end

