function [ data_filt ] = tukey_filt( data,coperiod,dt,costap_wid )
% JBR -- 10/21/16
% Apply cosine taper (Tukey window) in the frequency domain
%
% costap_wid: 0 => square window
%             1 => Hann window (looks similar to Gaussian)
%             Default is 0.5
%


% BUILD TUKEYWIN FILTER
Nt = length(data);
faxis = [0:(Nt-mod(Nt-1,2))/2 , -(Nt-mod(Nt,2))/2:-1]/dt/Nt;
fmax = 1/coperiod(1);
fmin = 1/coperiod(2);
[~ , If_pos] = find(faxis>= fmin & faxis<= fmax);
[~ , If_neg] = find(faxis<=-fmin & faxis>=-fmax);
costap_lenpos = length(If_pos);
costap_lenneg = length(If_neg);
costap_pos = tukeywin(costap_lenpos,costap_wid);
costap_neg = tukeywin(costap_lenneg,costap_wid);
costap_full = zeros(size(data));
costap_full(If_pos) = costap_pos;
costap_full(If_neg) = costap_neg;

% APPLY FILTER
data_filt = data.*costap_full;

if 0
    figure(99); clf;
    plot(1./faxis(1:Nt-1),abs(data(1:Nt-1))','-k'); hold on;
    plot(1./faxis(1:Nt-1),abs(data_filt(1:Nt-1))','-r');
    xlabel('Period');
    pause;
end


end

