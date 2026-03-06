% Example of windowing the fundamental mode in the Radon domain before
% transforming back to the time domain. In this example, we already know
% the true fundamental mode dispersion, so we build the windowing function
% around that. In practice, you will want to design your own windowing
% function around the energy you are interested in. In the future, we will
% add an interactive script for defining the window.
% 
% J. Russell - 3/2026
% github.com/jbrussell

clear;
setup_parameters;

% Bandpass filter for plotting purposes.
% There will likely be some issues at the edges of the Radon window, so filter to try to avoid those effects
f_min_filt = 1/120; % 1/150;
f_max_filt = 1/30; % 1/20;

% Load precalculated LRT
load([LRTmatpath,'LRT_',method,'_aug.mat']);
f = mat.f;
P_axis = mat.P_axis;
L = mat.L;
Rfft = mat.Rfft;
m_index = mat.m_index;

% Load PA5 dispersion
temp = load('./pa5_5km/dispersion_pa5_5km_b5.mat');
dat = temp.dat;
% Organize dispersion
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

% Apply bandpass filter for plotting waveforms 
M_filt = zeros(size(M));  
dt = t(2) - t(1);
for ii = 1:size(M_filt,1)
    dat_taper = cos_taper(M(ii,:)); 
    fs = 1/dt;
    [b,a] = butter(2,[f_min_filt/(fs/2) f_max_filt/(fs/2)]); % (20 - 150 seconds)
    %fvtool(b,a);
    M_filt(ii,:) = filtfilt(b,a,dat_taper);
end

%% Define Gaussian window function (FREQUENCY)

% Window around fundamental mode in the period band of interest
sigma_eff = .2; % control width of gaussian window
per_picks = DISP(1).Tq;
phv_picks = DISP(1).cvq;
per_v = 1./f;
phv_v = 1./P_axis;
per_v_win = per_v(per_v>=min(mat.per_vec) & per_v<=max(mat.per_vec));
phv_v_ref = interp1(per_picks,phv_picks,per_v_win);
phv_v_ref(isnan(phv_v_ref)) = 0;
sigma = sigma_eff*per_v_win/max(per_v_win);
for iper = 1:length(per_v_win)
    gauss(:,iper) = 1/(sigma(iper)*sqrt(2*pi))*exp(-1/2*((phv_v(:)-phv_v_ref(iper))/sigma(iper)).^2);
    gauss_norm(:,iper) = gauss(:,iper)/max(gauss(:,iper));
end
% Build augmented window matrix containing all frequencies (will be zero outside of window defined above)
gauss_aug = zeros(size(mat.Rfft));
I_posf = find(f>=min(1./per_v_win) & f<=max(1./per_v_win));
I_negf = find(f>=-1*max(1./per_v_win) & f<=-1*min(1./per_v_win));
gauss_aug(:,I_posf) = gauss_norm;
gauss_aug(:,I_negf) = flip(gauss_norm,2);


% Plot window function and 2-standard deviations
figure(13);
set(gcf,'Position',[56         291        1082         403]);
colormap([ones(30,3).*[0.2665 0.0033 0.3273]; viridis(100)]);
FS = 15;
subplot(1,2,1); box on; hold on;
h = pcolor(per_v_win,phv_v,gauss_norm); 
shading interp;
% set(h, 'edgecolor', 'none');
plot(per_v_win,phv_v_ref,'--','color',[.5 .5 .5],'linewidth',1.5);
plot(per_v_win,phv_v_ref+2*sigma,'--r','linewidth',1.5);
plot(per_v_win,phv_v_ref-2*sigma,'--r','linewidth',1.5);
colorbar;
set(gca,'YDir','normal','FontSize',FS,'linewidth',1.5,'TickDir','out');
title('Gaussian Window Function (Fund)');
xlabel('Period (s)');
ylabel('Phase Velocity (km/s)');
ylim([mat.v_min mat.v_max]);
xlim([min(mat.per_vec) max(mat.per_vec)]);

subplot(1,2,2); 
if is_globnorm
    imagesc(mat.per_vec, mat.phv_vec,  abs(mat.R_Tv)./prctile(mat.R_Tv(:),99)); hold on;
else
    imagesc(mat.per_vec, mat.phv_vec,  abs(mat.R_Tv)./max(abs(mat.R_Tv))); hold on;
end
for ii = 1:BRANCHES
    plot(DISP(ii).Tq(1:10:end),DISP(ii).cvq(1:10:end),'-','color',[1 0 0],'linewidth',1.5);   
end
plot(per_v_win,phv_v_ref+2*sigma,'--r','linewidth',1.5);
plot(per_v_win,phv_v_ref-2*sigma,'--r','linewidth',1.5);
caxis([0 1]);
xlim([min(mat.per_vec) max(mat.per_vec)]);
ylim([mat.v_min mat.v_max]);
title(method,'Interpreter','none'); ylabel('Velocity (km/s)'); xlabel('Period (s)');
set(gca,'YDir','normal','FontSize',FS,'linewidth',1.5,'TickDir','out');

%% Apply window in Radon domain and transform back to time (x-f)

% Define a sparse augmented Radon panel
Rfft_win = Rfft.*gauss_aug; 
Rfft_win_aug_sparse = repmat(Rfft_win,length(f),1).*m_index;
 
% Forward calculate data in time domain
M_win_fft = L*Rfft_win_aug_sparse; % d_est
M_win = ifft(M_win_fft,[],2); % Transform truncated original spectra back to time
M_win = real(M_win);

% Apply bandpass filter to windowed data 
M_win_filt = zeros(size(M_win));  
for ii = 1:size(M_win,1)
    dat_taper = cos_taper(M_win(ii,:)); 
    fs = 1/dt;
    [b,a] = butter(2,[f_min_filt/(fs/2) f_max_filt/(fs/2)]); % (20 - 150 seconds)
    %fvtool(b,a);
    M_win_filt(ii,:) = filtfilt(b,a,dat_taper);
end

%% Load true fundamental mode synthetics and compare
% Load fundamental mode waveforms
fundmode = load('./pa5_5km/Synth_120W_150W_fund.mat');
Delta_fund = deg2km(fundmode.Delta)';
% Apply bandpass filter for plotting waveforms 
M_filt_fund = zeros(size(fundmode.M));  
dt_fund = fundmode.t(2) - fundmode.t(1);
for ii = 1:size(M_filt_fund,1)
    dat_taper = cos_taper(fundmode.M(ii,:)); 
    fs = 1/dt_fund;
    [b,a] = butter(2,[f_min_filt/(fs/2) f_max_filt/(fs/2)]); % (20 - 150 seconds)
    %fvtool(b,a);
    M_filt_fund(ii,:) = filtfilt(b,a,dat_taper);
end

figure(14); clf;
amp = 30;
set(gcf,'position',[428          14        1025        1004],'color','w');
subplot(1,2,1); box on; hold on;
plot(t,M_filt./max(M_filt,[],2)*amp+Delta','-k','linewidth',1);
plot(t,M_filt_fund./max(M_filt_fund,[],2)*amp+Delta_fund','-b','linewidth',1);
xlabel('Time (s)'); ylabel('Distance (km)');
set(gca,'YDir','reverse','FontSize',FS,'linewidth',1.5);
xlim([400 1300]);
title(['Full mode data (',num2str(1./f_max_filt),'-',num2str(1./f_min_filt),' s)']);
h1(1) = plot(t,M_filt(1,:)/max(M_filt(1,:))*amp+Delta(1),'-k','linewidth',1);
h1(2) = plot(t,M_filt_fund(1,:)/max(M_filt_fund(1,:))*amp+Delta_fund(1),'-b','linewidth',1);
legend(h1,{'Data';'True fund only'},'location','northeast');

subplot(1,2,2); box on; hold on;
plot(t,M_filt_fund./max(M_filt_fund,[],2)*amp+Delta_fund','-b','linewidth',1);
plot(t,M_win_filt./max(M_win_filt,[],2)*amp+Delta','-r','linewidth',1);
xlabel('Time (s)'); ylabel('Distance (km)');
set(gca,'YDir','reverse','FontSize',FS,'linewidth',1.5);
xlim([400 1300]);
title(['Windowed data (',num2str(1./f_max_filt),'-',num2str(1./f_min_filt),' s)']);
h2(1) = plot(t,M_filt_fund(1,:)/max(M_filt_fund(1,:))*amp+Delta_fund(1),'-b','linewidth',1);
h2(2) = plot(t,M_win_filt(1,:)/max(M_win_filt(1,:))*amp+Delta(1),'-r','linewidth',1);
legend(h2,{'True fund only';'Windowed data'},'location','northeast');

% if ~exist(figpath)
%     mkdir(figpath);
% end
% save2pdf([figpath,'LRT_',method,'_isolate_modes.pdf'],3,500);
