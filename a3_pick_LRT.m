% Interactively pick dispersion from radon panel
% 
% J. Russell
% github.com/jbrussell

clear;
setup_parameters;

% Load precalculated LRT
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
    absR_Tv = abs(mat.R_Tv)./prctile(mat.R_Tv(:),99);
else
    absR_Tv = abs(mat.R_Tv)./max(abs(mat.R_Tv));
end
phv_trace = [];
per_trace = [];
phv_trace_std = [];
ipk = 0;
for iper = 1:Npers
    [~,I_per] = min(abs(mat.per_vec-pers(iper)));
    [pks,locs,w,p] = findpeaks(absR_Tv(:,I_per),mat.phv_vec,'MinPeakProminence',min_peak_prom,'MinPeakDistance',min_peak_dist);
    
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

%% MANUAL PICKING OF DISPERSION CURVE

[picks_LRT] = pick_LRT_disp(absR_Tv, mat, per_trace, phv_trace, phv_trace_std);

if ~exist(picks_out_path)
    mkdir(picks_out_path);
end
% Save picks to mat file
if ~isempty(picks_LRT)
    save([picks_out_path,'LRTpicks_',method,'.mat'],'picks_LRT');
end

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
    plot(DISP(ii).Tq(1:10:end),DISP(ii).cvq(1:10:end),'-','color',[0 0 0],'linewidth',1.5);   
end
caxis([0 1]);
xlim([min(mat.per_vec) max(mat.per_vec)]);
ylim([mat.v_min mat.v_max]);
title(method,'Interpreter','none'); ylabel('Velocity (km/s)'); xlabel('Period (s)');
set(gca,'YDir','normal','FontSize',FS,'linewidth',1.5,'TickDir','out');

colormap([ones(30,3).*[0.2665 0.0033 0.3273]; viridis(100)]);
pos = get(gca,'Position');
cb = colorbar;
set(cb,'linewidth',1.5);
set(gca,'Position',pos);
