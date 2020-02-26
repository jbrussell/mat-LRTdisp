% Plot the Radon panel calculated in a1
% 
% J. Russell
% github.com/jbrussell

clear;
setup_parameters;

% Load precalculated LRT
load([LRTmatpath,'LRT_',method,'_',comp,'.mat']);

% Load waveforms
load(mat.ndata);
Delta = Delta';

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
% [Delta_srt,I_srt] = sort(Delta);
% imagesc(t,Delta_srt,M(I_srt,:)); hold on;
plot(t,M./max(M,[],2)*10+Delta','-k','linewidth',1);
% title('Love waves (0T-4T)'); 
xlabel('Time (s)'); ylabel('Distance (km)');
set(gca,'YDir','reverse','FontSize',FS,'linewidth',1.5);

subplot(1,2,2); 
if is_globnorm
    imagesc(mat.per_vec, mat.phv_vec,  abs(mat.R_Tv)./prctile(mat.R_Tv(:),99)); hold on;
% contour(mat.per_vec, mat.phv_vec,  abs(mat.R_Tv)./prctile(mat.R_Tv(:),99)>=0.9,[1],'-w','linewidth',2); hold on;
else
    imagesc(mat.per_vec, mat.phv_vec,  abs(mat.R_Tv)./max(abs(mat.R_Tv))); hold on;
end
errorbar(per_trace,phv_trace,phv_trace_std,'or','MarkerFaceColor',[1 1 1],'linewidth',1.5,'markersize',7);
for ii = 1:BRANCHES
    plot(DISP(ii).Tq(1:10:end),DISP(ii).cvq(1:10:end),'-','color',[1 1 1],'linewidth',1.5);
    plot(DISP(ii).Tq(1:10:end),DISP(ii).cvq(1:10:end),'--','color',[0 0 0],'linewidth',1.5);
end
% colorbar;
caxis([0 1]);
% caxis([0 2e-3])
xlim([min(mat.per_vec) max(mat.per_vec)]);
% xlim([3 14]);
ylim([mat.v_min mat.v_max]);
title(method,'Interpreter','none'); ylabel('Velocity (km/s)'); xlabel('Period (s)');
set(gca,'YDir','normal','FontSize',FS,'linewidth',1.5,'TickDir','out');


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
% colormap([zeros(20,3); magma; ones(1,3)]);

pos = get(gca,'Position');
cb = colorbar;
set(cb,'linewidth',1.5);
set(gca,'Position',pos);
%%%%

if ~exist(figpath)
    mkdir(figpath);
end
save2pdf([figpath,'LRT_',method,'_',comp,'.pdf'],3,100);
