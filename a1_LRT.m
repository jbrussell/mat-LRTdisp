% Calculate radon transform
%
% t        - time axis.
% Delta    - distance (offset) axis.
% M        - Amplitudes of phase arrivals.
% indicies - list of indicies relevent to the S670S phase.
% 
% J. Russell
% github.com/jbrussell

clear;
setup_parameters;

% Save output?
is_savemat = 1;

% Load synthetic data
load(ndata,'-mat');
Delta = deg2km(Delta');
delta=mean(Delta);

% Load PA5 dispersion
load('./pa5_5km/dispersion_pa5_5km_b5.mat');

% Organize dipsersion
BRANCHES=5;
for ii = 1:BRANCHES
    DISP(ii).n = ii-1;
    DISP(ii).cv =  dat{ii}(:,6);
    DISP(ii).gv =  dat{ii}(:,7);
    DISP(ii).cvq = dat{ii}(:,8);
    DISP(ii).Tq =  dat{ii}(:,9);
    DISP(ii).T =   dat{ii}(:,10);
%     plot(Tq(1:10:end),cvq(1:10:end),'--','color',[.5 .5 .5],'linewidth',1);   
end

%% Do Radon Transform
% Invert to Radon domain using several different methods with varying
% degrees of sparseness
tic;
% [ Rfft,f ] = Radon_conjgrad(P_axis,t,M,Delta,maxiter,rthresh,method);
[ Rfft,f ] = Radon_conjgrad_fast(P_axis,t,M,Delta,f_min,f_max,maxiter,rthresh,method);
toc

[~,I_fmin_plot] = min(abs(f-f_min)); [~,I_fmax_plot] = min(abs(f-f_max));
I_fmin_plot = max(1, I_fmin_plot-1); I_fmax_plot = min(length(f), I_fmax_plot+1);
fplot = f(I_fmin_plot:I_fmax_plot);
[~,I_pmin_plot] = min(abs(P_axis-1/v_max)); [~,I_pmax_plot] = min(abs(P_axis-1/v_min)); 
I_pmin_plot = max(1, I_pmin_plot-1); I_pmax_plot = min(length(P_axis), I_pmax_plot+1);
P_axisplot = P_axis(I_pmin_plot:I_pmax_plot);

Rfftplot = Rfft(I_pmin_plot:I_pmax_plot,I_fmin_plot:I_fmax_plot);
% [ perplot,vplot,R_Tv ] = FreqSlow2PeriodVeloc( fplot,P_axisplot,abs(Rfftplot));
[ perplot,vplot,R_Tv ] = FreqSlow2PeriodVeloc( fplot,P_axisplot,Rfftplot);

% Convert from frequency-slowness to frequency-velocity
[Fplot, Vplot2, R_Fv] = FreqSlow2FreqVeloc(fplot, P_axisplot, Rfftplot);

% Apply bandpass filter for plotting waveforms 
dat_filt = zeros(size(M));  
dt = t(2) - t(1);
costap_wid = 0.2; % 0 => box filter; 1 => Hann window
for ii = 1:size(dat_filt,1)
    dat_taper = cos_taper(M(ii,:)); 
    [dat_filt_fft] = tukey_filt( fft(fftshift(dat_taper)),[1/f_max 1/f_min],dt,costap_wid);
    dat_filt(ii,:) = fftshift(real(ifft(dat_filt_fft)));
end

%%
% Plot figures.
figure(3); clf;
set(gcf,'Position',[173.0000  262.0000  880.0000  438.0000]);


subplot(1,2,1); hold on;
% plot(t,M./max(M,[],2)*50+Delta','-k','linewidth',1);
plot(t,dat_filt./max(dat_filt,[],2)*1+Delta','-k','linewidth',1);
% title('Love waves (0T-4T)'); 
xlabel('Time (s)'); ylabel('Distance (km)');
set(gca,'YDir','reverse');
xlim([400 1300]);

subplot(1,2,2); 
if is_globnorm
    imagesc(perplot(1,1:end), vplot(1:end,1),  abs(R_Tv)./prctile(abs(R_Tv(:)),99)); hold on;
% contour(mat.per_vec, mat.phv_vec,  abs(mat.R_Tv)./prctile(mat.R_Tv(:),99)>=0.9,[1],'-w','linewidth',2); hold on;
else
    imagesc(perplot(1,1:end), vplot(1:end,1),  abs(R_Tv)./max(abs(R_Tv))); hold on;
end
for ii = 1:BRANCHES
    plot(DISP(ii).Tq(1:10:end),DISP(ii).cvq(1:10:end),'-','color',[1 0 0],'linewidth',1.5);   
end
caxis([0 1]);
xlim([min(perplot(1,1:end)) max(perplot(1,1:end))]);
ylim([v_min v_max]);
title(method); ylabel('Velocity (km/s)'); xlabel('Period (s)');
set(gca,'YDir','normal');


%%%%%%%%%%%% make colormap %%%%%%%%%%%%
colormap([ones(30,3).*[0.2665 0.0033 0.3273]; viridis(100)]);
pos = get(gca,'Position');
cb = colorbar;
set(cb,'linewidth',1.5);
set(gca,'Position',pos);
%%%%

% figpath = './figs/';
% if ~exist(figpath)
%     mkdir(figpath);
% end
% % save2pdf([figpath,'LRT_',method,'_',comp,'.pdf'],3,100);

%% Save results to mat
if is_savemat
    if ~exist(LRTmatpath)
        mkdir(LRTmatpath);
    end
    mat.R_Tv = abs(R_Tv); % abs() to save space
    mat.per_vec = perplot(1,1:end);
    mat.phv_vec = vplot(1:end,1);
    
    % Fixed: Define frequency-velocity variables
    mat.R_Fv = abs(R_Fv);
    mat.freq_vec = Fplot(1,1:end);
    mat.phv_freq_vec = Vplot2(1:end,1);
    
    mat.ndata = ndata;
    mat.inversion.method = method;
    mat.inversion.maxiter = maxiter;
    mat.inversion.rthresh = rthresh;
    mat.f_min = f_min;
    mat.f_max = f_max;
    mat.v_min = v_min;
    mat.v_max = v_max;
    save([LRTmatpath,'LRT_',method,'.mat'],'mat');
end

%% Plot f-v
% Plot figures.
figure(4); clf;
set(gcf,'Position',[173.0000  262.0000  880.0000  438.0000]);


subplot(1,2,1); hold on;
% plot(t,dat./max(dat,[],2)*10+Delta','-k','linewidth',1);
plot(t,dat_filt./max(dat_filt,[],2)*1+Delta','-k','linewidth',1);
% title('Love waves (0T-4T)'); 
xlabel('Time (s)'); ylabel('Distance (km)');
set(gca,'YDir','reverse');

subplot(1,2,2); 
if is_globnorm
    imagesc(mat.freq_vec(1,1:end), mat.phv_freq_vec(1:end,1),  abs(mat.R_Fv)./prctile(abs(mat.R_Fv(:)),99)); hold on;
% contour(mat.per_vec, mat.phv_vec,  abs(mat.R_Tv)./prctile(mat.R_Tv(:),99)>=0.9,[1],'-w','linewidth',2); hold on;
else
    imagesc(mat.freq_vec(1,1:end), mat.phv_freq_vec(1:end,1),  abs(mat.R_Fv)./max(abs(mat.R_Fv))); hold on;
end

% colorbar;
caxis([0 1]);
% caxis([0 2e-3])
xlim([min(mat.freq_vec(1,1:end)) max(mat.freq_vec(1,1:end))]);
% xlim([3 14]);
ylim([v_min v_max]);
title(method); ylabel('Velocity (km/s)'); xlabel('Frequency (Hz)');
set(gca,'YDir','normal');
colormap([ones(30,3).*[0.2665 0.0033 0.3273]; viridis(100)]);
title(method); ylabel('Velocity (km/s)'); xlabel('Frequency (Hz)');
set(gca,'YDir','normal');
colormap([ones(30,3).*[0.2665 0.0033 0.3273]; viridis(100)]);
