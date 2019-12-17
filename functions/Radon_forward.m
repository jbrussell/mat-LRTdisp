function M=Radon_forward(t,p,R,delta,ref_dist,line_model)
%This function applies the time-shift Radon operator A, to the Radon 
%domain.  Will calculate the move-out data, given the inputs:
% -t        -- vector of time axis.
% -p        -- vector of slowness axis you would like to invert to.
% -R        -- matrix of Radon data, ordered size(R)==[length(p),length(t)].
% -delta    -- vector of distance axis.
% -ref_dist -- reference distance the path-function will shift about.
%
% -line_model, select one of the following options for path integration:
%     'linear'     - linear paths in the spatial domain (default)
%     'parabolic'  - parabolic paths in the spatial domain.
%
%Output spatial domain is ordered size(M)==[length(delta),length(t)].
%
%Known limitations:
% - Assumes evenly sampled time axis.
% - Assumes Radon data isn't complex.
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
  
  % Exit if inconsistent data is input.
  if(min([ip,it]~=size(R)))
      fprintf('Dimensions inconsistent!\nsize(R)~=[length(p),length(t)]\n');
      M=0;
      return;
  end;
  
  % Preallocate space in memory.
  Mfft=zeros(iDelta, iF);
  A=zeros(iDelta,ip);
  Tshift=A;

  % Define some values.
  Dist_array=delta-ref_dist;
  dF=1/(t(1)-t(2));
  Rfft=fft(R,iF,2);
  
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
      
      % Apply Radon operator.
      Mfft(:,i)=A*Rfft(:,i);
      
      % Assuming Hermitian symmetry of the fft make negative frequencies the complex conjugate of current solution.
      if(i ~= 1)
          Mfft(:,iF-i+2)=conj(Mfft(:,i));
      end
  end;
  
  M=ifft(Mfft,iF,2, 'symmetric');
  M=M(:,1:it);
return;
