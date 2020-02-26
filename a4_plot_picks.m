% Plot dispersion picks from a3
% 
% J. Russell
% github.com/jbrussell

clear;
setup_parameters;

% Load precalculated LRT
load([LRTmatpath,'LRT_',method,'_',comp,'.mat']);

% Load dispersion picks
load([picks_out_path,'LRTpicks_',method,'_',comp,'.mat']);

%%
% Plot picks  to double check

figure(4); clf;
FS = 15;
if is_globnorm
    imagesc(mat.per_vec, mat.phv_vec,  abs(mat.R_Tv)./prctile(mat.R_Tv(:),99)); hold on;
% contour(mat.per_vec, mat.phv_vec,  abs(mat.R_Tv)./prctile(mat.R_Tv(:),99)>=0.9,[1],'-w','linewidth',2); hold on;
else
    imagesc(mat.per_vec, mat.phv_vec,  abs(mat.R_Tv)./max(abs(mat.R_Tv))); hold on;
end
clr = lines(length(picks_LRT));
for itr = 1:length(picks_LRT)
    errorbar(picks_LRT(itr).per,picks_LRT(itr).phv,picks_LRT(itr).phv_std,'-o','Color',clr(itr,:),'LineWidth',1.5);
end
% errorbar(per_trace,phv_trace,phv_trace_std,'or','MarkerFaceColor',[1 1 1],'linewidth',1.5,'markersize',7);
for ii = 1:BRANCHES
    plot(DISP(ii).Tq(1:10:end),DISP(ii).cvq(1:10:end),'-','color',[1 0 0],'linewidth',1.5);
%     plot(DISP(ii).Tq(1:10:end),DISP(ii).cvq(1:10:end),'--','color',[0 0 0],'linewidth',1.5);
end
% colorbar;
caxis([0 1]);
% caxis([0 2e-3])
xlim([min(mat.per_vec) max(mat.per_vec)]);
% xlim([3 14]);
ylim([mat.v_min mat.v_max]);
title(method,'Interpreter','none'); ylabel('Velocity (km/s)'); xlabel('Period (s)');
set(gca,'YDir','normal','FontSize',FS,'linewidth',1.5,'TickDir','out');

colormap([ones(30,3).*[0.2665 0.0033 0.3273]; viridis(100)]);

pos = get(gca,'Position');
cb = colorbar;
set(cb,'linewidth',1.5);
set(gca,'Position',pos);
%%%%

save2pdf([figpath,'LRT_picks_',method,'_',comp,'.pdf'],4,100);
