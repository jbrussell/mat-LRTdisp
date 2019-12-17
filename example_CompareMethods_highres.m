% Compare Radon transform methods applied to synthetic multimode (0T-4T) 
% Love waves. The velocity model used to calculate the synethic waveforms in this
% example is PA5 (Gaherty et al. 1996). The synthetic station configuration
% is a linear array extending 30 degrees.
%
% So-called "High resolution linear radon transform" methods:
% - L1 (Schultz, 2012)
% - Cauchy (Schultz, 2012)
% - CG_IRLS (Iterative Reweighted Least-Squares; Li 2006)
% 
% J. Russell
% github.com/jbrussell

clear;

addpath('./functions/CG_methods/');

% Load variables.
load('./pa5_5km/Synth_120W_150W.mat','-mat');

isnoise = 1; % add gaussian noise?

Delta = deg2km(Delta');
% t        - time axis.
% Delta    - distance (offset) axis.
% M        - Amplitudes of phase arrivals.
% indicies - list of indicies relevent to the S670S phase.

load('./pa5_5km/dispersion_pa5_5km_b5.mat');
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

% Define some variables for RT.
maxiter = 10; %100;
rthresh = 1e-6;
method1='L2'; mu1=5e-2; %5e-2; % L2
method2='L1'; mu2=[5e-2 5e-6]; % L1
method3='Cauchy'; mu3=[5e-2 1e-8]; % Cauchy
P_axis=10:0.1:40; % s/deg
P_axis = P_axis / 111; %(s/km);
delta=mean(Delta);
f_min = 1/150;
f_max = 1/20;
v_min = 4;
v_max = 8;

% Add noise to data
std_noise = rms(M(:))*0.5;
if isnoise
    M = M + normrnd(0,std_noise,size(M));
end

% Invert to Radon domain using several different methods with varying
% degrees of sparseness
tic;
[R,Rfft,f]=Radon_inverse_disp(t, Delta, M, P_axis, ones(size(Delta)), delta, 'Linear', 'L2', mu1);
toc

tic;
[~,Rfft_2,f_2 ] = Radon_inverse_disp(t, Delta, M, P_axis, ones(size(Delta)), delta, 'Linear', 'L1', mu2);
toc

tic;
[~,Rfft_3,f_3 ] = Radon_inverse_disp(t, Delta, M, P_axis, ones(size(Delta)), delta, 'Linear', 'Cauchy', mu3);
toc

tic;
[ Rfft_4,f_4 ] = Radon_conjgrad(P_axis,t,M,Delta,maxiter,rthresh,'CG_IRLS');
toc

[~,I_fmin_plot] = min(abs(f-f_min)); [~,I_fmax_plot] = min(abs(f-f_max));
I_fmin_plot=I_fmin_plot-1; I_fmax_plot=I_fmax_plot+1;
fplot = f(I_fmin_plot:I_fmax_plot);
[~,I_pmin_plot] = min(abs(P_axis-1/v_max)); [~,I_pmax_plot] = min(abs(P_axis-1/v_min)); 
I_pmin_plot=I_pmin_plot-1; I_pmax_plot=I_pmax_plot+1;
P_axisplot = P_axis(I_pmin_plot:I_pmax_plot);

Rfftplot = Rfft(I_pmin_plot:I_pmax_plot,I_fmin_plot:I_fmax_plot);
[ perplot,vplot,R_Tv ] = FreqSlow2PeriodVeloc( fplot,P_axisplot,abs(Rfftplot));
[ ~,~,R_Tv_2 ] = FreqSlow2PeriodVeloc( fplot,P_axisplot,abs(Rfft_2(I_pmin_plot:I_pmax_plot,I_fmin_plot:I_fmax_plot)));
[ ~,~,R_Tv_3 ] = FreqSlow2PeriodVeloc( fplot,P_axisplot,abs(Rfft_3(I_pmin_plot:I_pmax_plot,I_fmin_plot:I_fmax_plot)));
[ ~,~,R_Tv_4 ] = FreqSlow2PeriodVeloc( fplot,P_axisplot,abs(Rfft_4(I_pmin_plot:I_pmax_plot,I_fmin_plot:I_fmax_plot)));


%%
% Plot figures.
figure(3); clf;


subplot(2,3,1);
[Delta_srt,I_srt] = sort(Delta);
imagesc(t,Delta_srt,M(I_srt,:)); hold on;
for ii = 1:length(Delta)
    plot(t,M(ii,:)/max(M(ii,:))*30+Delta(ii),'-k','linewidth',1);
end
title('Love waves (0T-4T)'); xlabel('Time (s)'); ylabel('Distance (km)');

subplot(232); imagesc(perplot(1,1:end), vplot(1:end,1),  abs(R_Tv(:,1:length(perplot))')'./max(abs(R_Tv(:,1:length(perplot))')')); hold on;
for ii = 1:BRANCHES
    plot(DISP(ii).Tq(1:10:end),DISP(ii).cvq(1:10:end),'-','color',[0 0.8 0],'linewidth',1);   
end
% colorbar;
caxis([0 1]);
xlim([min(perplot(1,1:end)) max(perplot(1,1:end))]);
ylim([v_min v_max]);
title('L2'); ylabel('Velocity (km/s)'); xlabel('Period (s)');
set(gca,'YDir','normal');

subplot(233); imagesc(perplot(1,1:end), vplot(1:end,1),  abs(R_Tv_2(:,1:length(perplot))')'./max(abs(R_Tv_2(:,1:length(perplot))')')); hold on;
for ii = 1:BRANCHES
    plot(DISP(ii).Tq(1:10:end),DISP(ii).cvq(1:10:end),'-','color',[0 0.8 0],'linewidth',1);   
end
% colorbar;
caxis([0 1]);
xlim([min(perplot(1,1:end)) max(perplot(1,1:end))]);
ylim([v_min v_max]);
title('L1'); ylabel('Velocity (km/s)'); xlabel('Period (s)');
set(gca,'YDir','normal');

subplot(234); imagesc(perplot(1,1:end), vplot(1:end,1),  abs(R_Tv_3(:,1:length(perplot))')'./max(abs(R_Tv_3(:,1:length(perplot))')')); hold on;
for ii = 1:BRANCHES
    plot(DISP(ii).Tq(1:10:end),DISP(ii).cvq(1:10:end),'-','color',[0 0.8 0],'linewidth',1);   
end
% colorbar;
caxis([0 1]);
xlim([min(perplot(1,1:end)) max(perplot(1,1:end))]);
ylim([v_min v_max]);
title('Cauchy'); ylabel('Velocity (km/s)'); xlabel('Period (s)');
set(gca,'YDir','normal');

subplot(235); imagesc(perplot(1,1:end), vplot(1:end,1),  abs(R_Tv_4(:,1:length(perplot))')'./max(abs(R_Tv_4(:,1:length(perplot))')')); hold on;
for ii = 1:BRANCHES
    plot(DISP(ii).Tq(1:10:end),DISP(ii).cvq(1:10:end),'-','color',[0 0.8 0],'linewidth',1);   
end
% colorbar;
caxis([0 1]);
xlim([min(perplot(1,1:end)) max(perplot(1,1:end))]);
ylim([v_min v_max]);
title('CG IRLS'); ylabel('Velocity (km/s)'); xlabel('Period (s)');
set(gca,'YDir','normal');

%%%%%%%%%%%% make colormap %%%%%%%%%%%%
colmap1 = cmap('steelblue',10,10,15);
colmap2 = cmap('orange',6,40,20);
colmap3 = rgb('orangered');
colmap4 = cmap('red',8,1,45);
colmap = [colmap1;flipud(colmap2);colmap3;flipud(colmap4)];
colormap(colmap);