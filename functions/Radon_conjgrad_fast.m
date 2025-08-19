function [ Rfft,f ] = Radon_conjgrad_fast(p,t,M,delta,f_min,f_max,maxiter,rthresh,method)
% Solve radon transform using the conjugate gradient method. This version only calculates 
% the frequencies of interest defined by [f_min, f_max], which can be significantly faster.
% The disadvantage of this approach is that the waveforms cannot be reconstructed without 
% the full frequency axis.
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

% Make sure waveforms are tapered
for ii = 1:size(M,1)
    M(ii,:) = cos_taper(M(ii,:));
end

% Define some array/matrices lengths.
it=length(t);
iF=pow2(nextpow2(it)+1); % Double length
iDelta=length(delta);
ip=length(p);
Mfft=fft(M,iF,2);
dF=1/(t(1)-t(2));
f = ((  [1:floor((iF+1)/2)]  -1)/iF)*dF*-1;

f_ind = find(f>=f_min*0.9 & f<=f_max*1.1);
iF_ind = length(f_ind)*2;

% Define blocks
blk_w = ip; %block width within L
blk_h = iDelta; % block height within L
delta_block = repmat(delta,ip,1)'; %repmat(delta,1,blk_w);
p_block = repmat(p,iDelta,1); %repmat(p',blk_h,1);

m0 = zeros(length(p),1);
Rfft = zeros(ip,iF_ind);
rfft = zeros(size(Mfft,1),iF_ind);
ii = 0;
for j = f_ind(:)'
    ii = ii + 1;
    if mod(ii,250) == 0 || ii==1
        disp([num2str(ii),'/',num2str(length(f_ind))]);
    end
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
    Rfft(:,ii) = m;
    rfft(:,ii) = r;
%     Rfft_est(:,j) = L*m; % d_est
end

f = f(f_ind);
    