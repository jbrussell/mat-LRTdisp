function [R,Rfft,fvec]=Radon_inverse_disp(t,delta,M,p,weights,ref_dist,line_model,inversion_model,hyperparameters)
% J. Russell
% Modified version to also output the Radon panel in the frequency
% domain, Rfft.
%
%This function inverts move-out data to the Radon domain given the inputs:
% -t        -- vector of time axis.
% -delta    -- vector of distance axis.
% -M        -- matrix of move-out data, ordered size(M)==[length(delta),length(t)].
% -p        -- vector of slowness axis you would like to invert to.
% -weights  -- weighting vector that determines importance of each trace.
%              set vector to ones for no preference.
% -ref_dist -- reference distance the path-function will shift about.
%
% -line_model, select one of the following options for path integration:
%     'linear'     - linear paths in the spatial domain (default)
%     'parabolic'  - parabolic paths in the spatial domain.
%
% -inversion model, select one of the following options for regularization schema:
%     'L2'       - Regularized on the L2 norm of the Radon domain (default)
%     'L1'       - Non-linear regularization based on L1 norm and iterative
%                  reweighted least sqaures (IRLS) see Sacchi 1997.
%     'Cauchy'   - Non-linear regularization see Sacchi & Ulrych 1995
%
% -hyperparameters, trades-off between fitting the data and chosen damping.
%
%Output radon domain is ordered size(R)==[length(p),length(t)].
%
%Known limitations:
% - Assumes evenly sampled time axis.
% - Assumes move-out data isn't complex.
%
%
% References: Schultz, R., Gu, Y. J., 2012. Flexible, inversion-based Matlab 
%             implementation of the Radon Transform.  Computers and 
%             Geosciences [In Preparation]
%
%             An, Y., Gu, Y. J., Sacchi, M., 2007. Imaging mantle 
%             discontinuities using least-squares Radon transform. 
%             Journal of Geophysical Research 112, B10303.
%
% Author: R. Schultz, 2012
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details: http://www.gnu.org/licenses/
%
  
  % Define some array/matrices lengths.
  it=length(t);
  iF=pow2(nextpow2(it)+1); % Double length
  iDelta=length(delta);
  ip=length(p);
  iw=length(weights);
  
  % Exit if inconsistent data is input.
  if(min([iDelta,it]~=size(M)))
      fprintf('Dimensions inconsistent!\nsize(M)~=[length(delta),length(t)]\n');
      R=0;
      return;
  end;
  if(iw~=iDelta)
      fprintf('Dimensions inconsistent!\nlength(delta)~=length(weights)\n');
      R=0;
      return;
  end;
  
  % Exit if improper hyperparameters are entered.
  if(strcmpi(inversion_model, 'L1')||strcmpi(inversion_model, 'Cauchy'))
      if(length(hyperparameters)~=2)
          fprintf('Improper number of trade-off parameters\n');
          R=0;
          return;
      end;
  else % The code's default is L2 inversion.
      if(length(hyperparameters)~=1)
          fprintf('Improper number of trade-off parameters\n');
          R=0;
          return;
      end;
  end;
  
  % Preallocate space in memory.
  R=zeros(ip, it); %#ok<NASGU>
  Rfft=zeros(ip, iF);
  A=zeros(iDelta,ip);
  Tshift=A;
  AtA=zeros(ip,ip); %#ok<NASGU>
  AtM=zeros(ip,1); %#ok<NASGU>
  Ident=speye(ip);

  % Define some values.
  Dist_array=delta-ref_dist;
  dF=1/(t(1)-t(2));
  Mfft=fft(M,iF,2);
  W=spdiags(weights',0,iDelta,iDelta);
  
  dCOST=0; %#ok<NASGU>
  COST_cur=0;
  COST_prev=0;
  
  % Populate ray parameter then distance data in time shift matrix.
  for j=1:iDelta
      if( strcmpi(line_model,'parabolic') )
          Tshift(j,:)=p;
      else % Linear is default.
          Tshift(j,:)=p;
      end;
  end;
  for k=1:ip
      if( strcmpi(line_model,'parabolic') )
          Tshift(:,k)=(2*ref_dist*Tshift(:,k).*Dist_array')+(Tshift(:,k).*(Dist_array.^2)');
      else % Linear is default.
          Tshift(:,k)=Tshift(:,k).*Dist_array';
      end;
  end;
  
  % Loop through each frequency.
  for i=1:floor((iF+1)/2)
      
      % Make time-shift matrix, A.
      f=((i-1)/iF)*dF;
      A=exp(  (2i*pi*f).*Tshift  );
      fvec(i) = -f;
      
      % M = A R  --->  AtM = AtA R
      % Solve the weighted, L2 least-squares problem for an initial solution.
      AtA=A'*W*A;
      AtM=A'*W*Mfft(:,i);
      mu=abs(trace(AtA))*hyperparameters(1);
      Rfft(:,i)=(AtA+mu*Ident)\AtM;
      
      % Non-quadratic inversions use IRLS to solve, iterate until solution convergence.
      if(strcmpi(inversion_model, 'Cauchy')||strcmpi(inversion_model, 'L1'))
          
          % Initialize hyperparameters.
          b=hyperparameters(2);
          lambda=mu*b;
          
          % Initialize cost functions.
          dCOST=inf;
          if(strcmpi(inversion_model, 'Cauchy'))
              COST_prev=norm(Mfft(:,i)-A*Rfft(:,i), 2)+lambda*sum(log( abs(Rfft(:,i) ).^2 + b ));
          elseif(strcmpi(inversion_model, 'L1'))
              COST_prev=norm(Mfft(:,i)-A*Rfft(:,i), 1)+lambda*norm(abs(Rfft(:,i))+b, 2);
          end;
          iter=1;
    
          % Iterate until negligible change to cost function.
          while( (dCOST > 0.001)&&(iter<20) )
              
              % Setup inverse problem.
              if(strcmpi(inversion_model, 'Cauchy'))
                  Q=spdiags( 1./( abs(Rfft(:,i) ).^2 + b), 0, ip, ip  );
              elseif(strcmpi(inversion_model, 'L1'))
                  Q=spdiags( 1./( abs(Rfft(:,i) ) + b), 0, ip, ip  );
              end;
              Rfft(:,i)=( lambda*Q+AtA )\AtM;
              
              % Determine change to cost function.
              if(strcmpi(inversion_model, 'Cauchy'))
                  COST_cur=norm(Mfft(:,i)-A*Rfft(:,i), 2)+lambda*sum(log( abs(Rfft(:,i) ).^2 + b )-log(b));
              elseif(strcmpi(inversion_model, 'L1'))
                  COST_cur=norm(Mfft(:,i)-A*Rfft(:,i), 1)+lambda*norm(abs(Rfft(:,i))+b, 2);
              end;
              dCOST=2*abs(COST_cur-COST_prev)/(abs(COST_cur)+abs(COST_prev));
              COST_prev=COST_cur;
              
              iter=iter+1;
          end;
      end;
      
      % Assuming Hermitian symmetry of the fft make negative frequencies the complex conjugate of current solution.
      if(i ~= 1)
          Rfft(:,iF-i+2)=conj(Rfft(:,i));
      end
  end;
  
  R=ifft(Rfft,iF,2, 'symmetric');
  R=R(:,1:it);
return;
