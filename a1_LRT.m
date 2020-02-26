% Compute the Linear Radon Transform (LRT)
% 
% J. Russell
% github.com/jbrussell

clear;
setup_parameters;

is_savemat = 1;

% Load data
load(ndata);
delta=mean(Delta);
Delta = Delta';

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
plot(t,M./max(M,[],2)*10+Delta','-k','linewidth',1);
% title('Love waves (0T-4T)'); 
xlabel('Time (s)'); ylabel('Distance (km)');
set(gca,'YDir','reverse');

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
% colorbar;
caxis([0 1]);
% caxis([0 2e-3])
xlim([min(perplot(1,1:end)) max(perplot(1,1:end))]);
% xlim([3 14]);
ylim([v_min v_max]);
title(method); ylabel('Velocity (km/s)'); xlabel('Period (s)');
set(gca,'YDir','normal');


%%%%%%%%%%%% make colormap %%%%%%%%%%%%
colmap1 = cmap('steelblue',10,10,15);
colmap2 = cmap('orange',6,40,20);
colmap3 = rgb('orangered');
colmap4 = cmap('red',8,1,45);
colmap = [colmap1;flipud(colmap2);colmap3;flipud(colmap4)];
% colormap(colmap);
% colormap(cptcmap('GMT_ocean'));
% colormap(cptcmap('GMT_haxby'));

% pos01 = [0    1    0.4 0.4001 0.5999   0.6000 1.0000];
% clrs = [193 0 0; 224 0 0; 252 214 0; 235 235 235; 235 235 235; 81 211 255; 0 73 209]/255;
% mycolormap = customcolormap(pos01, clrs, 200);
% colormap(mycolormap);

% colormap(radon_cmap(200));
colormap([ones(30,3).*[0.2665 0.0033 0.3273]; viridis(100)]);
% colorbar;
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
    mat.ndata = ndata;
    mat.inversion.method = method;
    mat.inversion.maxiter = maxiter;
    mat.inversion.rthresh = rthresh;
    mat.f_min = f_min;
    mat.f_max = f_max;
    mat.v_min = v_min;
    mat.v_max = v_max;
    save([LRTmatpath,'LRT_',method,'_',comp,'.mat'],'mat');
end