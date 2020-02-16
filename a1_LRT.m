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

addpath('./functions/'); addpath('./functions/CG_methods/');

% Load Synthetic Love wave data
ndata = './pa5_5km/Synth_120W_150W.mat';
load(ndata,'-mat');
Delta = deg2km(Delta');
% Load PA5 dispersion
load('./pa5_5km/dispersion_pa5_5km_b5.mat');

% Save output?
is_savemat = 1;

% Normalization option for plotting
is_globnorm = 1; % 1 for normalize radon panel by global max; 0 for column norm

% Define some variables for RT.
maxiter = 10; %100;
rthresh = 1e-6;
method = 'CGG_weight';
% method = 'CG_IRLS';
delta=mean(Delta);
f_min = 1/150;
f_max = 1/20;
v_min = 4;
v_max = 8;
P_axis = [111/(v_max*1.1) : 0.1 : 111/(v_min*0.9)]; % s/deg
P_axis = P_axis / 111; %(s/km);

isnoise = 0; % add gaussian noise?

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
[ Rfft,f ] = Radon_conjgrad(P_axis,t,M,Delta,maxiter,rthresh,method);
toc

[~,I_fmin_plot] = min(abs(f-f_min)); [~,I_fmax_plot] = min(abs(f-f_max));
I_fmin_plot=I_fmin_plot-1; I_fmax_plot=I_fmax_plot+1;
fplot = f(I_fmin_plot:I_fmax_plot);
[~,I_pmin_plot] = min(abs(P_axis-1/v_max)); [~,I_pmax_plot] = min(abs(P_axis-1/v_min)); 
I_pmin_plot=I_pmin_plot-1; I_pmax_plot=I_pmax_plot+1;
P_axisplot = P_axis(I_pmin_plot:I_pmax_plot);

Rfftplot = Rfft(I_pmin_plot:I_pmax_plot,I_fmin_plot:I_fmax_plot);
% [ perplot,vplot,R_Tv ] = FreqSlow2PeriodVeloc( fplot,P_axisplot,abs(Rfftplot));
[ perplot,vplot,R_Tv ] = FreqSlow2PeriodVeloc( fplot,P_axisplot,Rfftplot);


%%
% Plot figures.
figure(3); clf;
set(gcf,'Position',[173.0000  262.0000  880.0000  438.0000]);


subplot(1,2,1); hold on;
plot(t,M./max(M,[],2)*50+Delta','-k','linewidth',1);
title('Love waves (0T-4T)'); 
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
    LRTmatpath = './LRT_mats/';
    if ~exist(LRTmatpath)
        mkdir(LRTmatpath);
    end
    mat.R_Tv = abs(R_Tv); % abs() to save space
    mat.per_vec = perplot(1,1:end);
    mat.phv_vec = vplot(1:end,1);
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