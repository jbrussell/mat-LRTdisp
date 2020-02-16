% Plotting radon panel and test dispersion curve tracing parameters
% 
% J. Russell
% github.com/jbrussell

clear;

addpath('./functions/'); addpath('./functions/CG_methods/');

method = 'CGG_weight';
% method = 'CG_IRLS';

% Load precalculated LRT
LRTmatpath = './LRT_mats/';
load([LRTmatpath,'LRT_',method,'.mat']);

% Load PA5 dispersion
load('./pa5_5km/dispersion_pa5_5km_b5.mat');

% Normalization option for plotting
is_globnorm = 1; % 1 for normalize radon panel by global max; 0 for column norm

% Parameters for tracing dispersion curves
min_peak_prom = 0.3; % Minimum peak prominence, threshold for peak height
min_peak_dist = 0.1; % Minimum separation between chosen peaks [km/s]
Npers = 25; % Number for periods
pers = logspace(log10(20),log10(150),Npers); % period vector 

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

% Load waveforms
load(mat.ndata);
Delta = deg2km(Delta)';

%% Find peaks
if is_globnorm
    R_Tv = abs(mat.R_Tv)./prctile(mat.R_Tv(:),99);
else
    R_Tv = abs(mat.R_Tv)./max(abs(mat.R_Tv));
end
phv_trace = [];
per_trace = [];
phv_trace_std = [];
ipk = 0;
for iper = 1:Npers
    [~,I_per] = min(abs(mat.per_vec-pers(iper)));
    [pks,locs,w,p] = findpeaks(R_Tv(:,I_per),mat.phv_vec,'MinPeakProminence',min_peak_prom,'MinPeakDistance',min_peak_dist);
    
    for ii = 1:length(pks)
        ipk = ipk+1;
        phv_trace(ipk) = locs(ii);
        per_trace(ipk) = pers(iper);
        phv_trace_std(ipk) = w(ii)*0.2; % 20% of peak corresponds to 0.90
    end
    
    if 0
        figure(99); clf; hold on;
        findpeaks(R_Tv(:,I_per),mat.phv_vec,'MinPeakProminence',min_peak_prom,'MinPeakDistance',min_peak_dist,'Annotate','extents');
%         plot(mat.phv_vec,R_Tv(:,I_per)); hold on;
        plot(locs,pks,'or');
        pause;
    end
end

%%
% Plot figures.
figure(3); clf;
set(gcf,'Position',[54         292        1069         405]);
FS = 15;

subplot(1,2,1); box on; hold on;
plot(t,M./max(M,[],2)*200+Delta','-k','linewidth',1);
title('Love waves (0T-4T)'); 
xlabel('Time (s)'); ylabel('Distance (km)');
set(gca,'YDir','reverse','FontSize',FS,'linewidth',1.5);
xlim([400 1300]);

subplot(1,2,2); 
if is_globnorm
    imagesc(mat.per_vec, mat.phv_vec,  abs(mat.R_Tv)./prctile(mat.R_Tv(:),99)); hold on;
% contour(mat.per_vec, mat.phv_vec,  abs(mat.R_Tv)./prctile(mat.R_Tv(:),99)>=0.9,[1],'-w','linewidth',2); hold on;
else
    imagesc(mat.per_vec, mat.phv_vec,  abs(mat.R_Tv)./max(abs(mat.R_Tv))); hold on;
end
% errorbar(per_trace,phv_trace,phv_trace_std,'or','MarkerFaceColor',[1 1 1],'linewidth',1.5,'markersize',7);
for ii = 1:BRANCHES
    plot(DISP(ii).Tq(1:10:end),DISP(ii).cvq(1:10:end),'-','color',[1 0 0],'linewidth',1.5);   
end
caxis([0 1]);
xlim([min(mat.per_vec) max(mat.per_vec)]);
ylim([mat.v_min mat.v_max]);
title(method,'Interpreter','none'); ylabel('Velocity (km/s)'); xlabel('Period (s)');
set(gca,'YDir','normal','FontSize',FS,'linewidth',1.5,'TickDir','out');


%%%%%%%%%%%% make colormap %%%%%%%%%%%%
colormap([ones(30,3).*[0.2665 0.0033 0.3273]; viridis(100)]);

pos = get(gca,'Position');
cb = colorbar;
set(cb,'linewidth',1.5);
set(gca,'Position',pos);
%%%%

figpath = './figs/';
if ~exist(figpath)
    mkdir(figpath);
end
save2pdf([figpath,'LRT_',method,'.pdf'],3,500);
