% Plot the cross-spectra in the time domain for the individual station pairs. 
% Filter is built using a Tukey taper with sharpness controlled by costap_wid. 
% The typical butterworth filter is not precise enough for the higest frequencies
%
% https://github.com/jbrussell

clear;
setup_parameters_MATnoise;
IsFigure = 0;
IsFigure_GAUS = 0; % Plot frequency domain filtered and unfiltered

%======================= PARAMETERS =======================%
comp = 'ZZ'; %'ZZ'; %'RR'; %'TT';
coperiod = [3 10]; % Periods to filter between
amp = 8e0;
windir = 'window3hr';
windir_for_SNR = 'window3hr'; % Data to use for calculating SNR threshold (for plotting purposes)
trace_space = 0; % km
snr_thresh = 2;
dep_tol = [0 0]; % [sta1, sta2] OBS Depth tolerance;
max_grv = inf; %5.5;
min_grv = 0.7; %1.6
xlims = [-500 500];
ylims = [0 450];
IsButterworth = 1;

ccf_path = '/data/irma6/jrussel/YoungPacificORCA/ccf_raw/';
% ccf_path = '/data/irma6/jrussel/YoungPacificORCA/ccf_FTN/';

%%% --- Parameters to build up gaussian filters --- %%% 
% (effects the width of the filter in the frequency domain)
costap_wid = 0.2; % 0 => box filter; 1 => Hann window

isplotwin = 1; %1;
isploth20 = 0;
isfigure_snr = 0;

h20_grv = 1.5;
%==========================================================%

dt = parameters.dt;
stalist = parameters.stalist;
nsta = parameters.nsta;
nsta = length(stalist);
winlength = parameters.winlength;
figpath = parameters.figpath;
%fig_winlength_path = [figpath,'window',num2str(winlength),'hr/fullStack/'];
% custom directory names
    fig_winlength_path = [figpath,windir,'/fullStack/'];

%------------ PATH INFORMATION -------------%
% ccf_path = parameters.ccfpath;
%ccf_winlength_path = [ccf_path,'window',num2str(winlength),'hr/'];
    ccf_winlength_path = [ccf_path,windir,'/'];
ccf_singlestack_path = [ccf_winlength_path,'single/'];
ccf_daystack_path = [ccf_winlength_path,'dayStack/'];
ccf_monthstack_path = [ccf_winlength_path,'monthStack/'];
ccf_fullstack_path = [ccf_winlength_path,'fullStack/'];

ccf_stack_path = ccf_fullstack_path;

figpath = [fig_winlength_path,num2str(coperiod(1)),'_',num2str(coperiod(2)),'s/'];
% create figure directory
if ~exist(fig_winlength_path)
    mkdir(fig_winlength_path)
end
if ~exist(figpath)
    mkdir(figpath)
end

%% Load Depths
STAS = stalist;
LATS = stalat;
LONS = stalon;
DEPTHS = staz;

%%
ccf_path = [ccf_stack_path,'ccf',comp,'/',];
ccf_path_SNR = [parameters.ccfpath,windir_for_SNR,'/','fullStack/','ccf',comp,'/',];
npairall = 0;
%------------ LOAD DATA AND PLOT IN TIME DOMAIN -------------%
for ista1=1:nsta % loop over all stations
    sta1=char(stalist(ista1,:));
    sta1dir=[ccf_path,sta1]; % dir to have all cross terms about this central station
    sta1dir_SNR = [ccf_path_SNR,sta1];
    nstapair = 0;
    for ista2 = 1: nsta % loop over station pairs
        sta2 = char(stalist(ista2,:));
        
        % if same station, skip
        if(strcmp(sta1,sta2))
            continue
        end
                
        filename = sprintf('%s/%s_%s_f.mat',sta1dir,sta1,sta2);
        filename_SNR = sprintf('%s/%s_%s_f.mat',sta1dir_SNR,sta1,sta2);
        
        if ~exist(filename,'file') % check that ccf file exists
            disp(['not exist ',filename])
            continue;
        end
        nstapair = nstapair + 1;
        
        % Want sta1 to be closest to the coast so waves at -lag travel
        % towards coast and waves at +lag travel away from coast.
        filename_sta1sta2 = filename;
        filename_sta2sta1 = [ccf_path,sta2,'/',sta2,'_',sta1,'_f.mat'];
        test = load(filename_sta1sta2);
        sta1lon = test.stapairsinfo.lons(1);
        sta2lon = test.stapairsinfo.lons(2);
        if sta1lon > sta2lon
            filename = filename_sta2sta1;
        else
            filename = filename_sta1sta2;
        end
        
        %----------- LOAD DATA -------------%
        data = load(filename);
        data_SNR = load(filename_SNR);
        ccf = data.coh_sum./data.coh_num;
        ccf(isnan(ccf)) = 0;
        ccf_SNR = data_SNR.coh_sum./data_SNR.coh_num;
        ccf_SNR(isnan(ccf_SNR)) = 0;
        if size(ccf,1) == 1
            ccf = ccf';
            ccf_SNR = ccf_SNR';
        end
        
        %%
        
        if IsFigure_GAUS
            T = length(ccf_filtered);
            faxis = [0:1/T:1/dt/2,-1/dt/2+1/T:1/T:-1/T];
            ind = find(faxis>0);
            figure(1); clf;
            plot((faxis(ind)),abs(ccf(ind)),'-k','linewidth',4); hold on;
            plot((faxis(ind)),abs(ccf_filt(ind)),'r');
            xlabel('Frequency (Hz)');
            title(['Gaussian filter ',sta1,'-',sta2,' (',num2str(coperiod(1)),'-',num2str(coperiod(2)),' s)']);
            pause;
        end
        %----------- Frequency ==> Time domain -------------%
        N = length(ccf);
        ccf_ifft = real(ifft(ccf,N)); % inverse FFT to get time domain
        ccf_ifft = fftshift(ccf_ifft); % rearrange values as [-lag lag]
        ccf_ifft = detrend(ccf_ifft);
        ccf_ifft = cos_taper(ccf_ifft);
        ccf_SNR_ifft = real(ifft(ccf_SNR,N)); % inverse FFT to get time domain
        ccf_SNR_ifft = fftshift(ccf_SNR_ifft); % rearrange values as [-lag lag]
        ccf_SNR_ifft = detrend(ccf_SNR_ifft);
        ccf_SNR_ifft = cos_taper(ccf_SNR_ifft);
        
        %----------- FILTER DATA (FREQUENCY DOMAIN) -------------%
        f1 = 1/coperiod(2);
        f2 = 1/coperiod(1);
        
        if ~IsButterworth
            [ ccf_filtered ] = tukey_filt( fft(fftshift(ccf_ifft)),coperiod,dt,costap_wid );
            [ ccf_filtered_SNR ] = tukey_filt(fft(fftshift(ccf_SNR_ifft)),coperiod,dt,costap_wid );
            ccf_ifft = fftshift(real(ifft(ccf_filtered)));
        else
            % Do butterworth filtering after rearranging
            [b, a] = butter(2,[f1 f2]*2*dt); % Butterworth Filter
            ccf_ifft = filtfilt(b,a,ccf_ifft);
            ccf_SNR_ifft = filtfilt(b,a,ccf_SNR_ifft);
            ccf_filtered = fft(fftshift(ccf_ifft));
            ccf_filtered_SNR = fft(fftshift(ccf_SNR_ifft));
            
            if 0
                figure(99); clf;
                [h,f] = freqz(b,a,length(ccf_ifft),dt);
                subplot(2,1,1);
                plot(f,abs(h));
                subplot(2,1,2);
                plot(f,angle(h)*180/pi);
            end
        end
        
        ccf_filt{nstapair} = ccf_ifft;
        
        %----------- NORMALIZE CCF FUNCTION -------------%
        ccf_filt{nstapair} = ccf_filt{nstapair}/max(abs(ccf_filt{nstapair}));
        
        % Distance between sta1 and sta2
        sta1sta2_dist(nstapair) = deg2km(distance(data.stapairsinfo.lats(1),data.stapairsinfo.lons(1),data.stapairsinfo.lats(2),data.stapairsinfo.lons(2)));
        stalats(ista1) = data.stapairsinfo.lats(1);
        stalons(ista1) = data.stapairsinfo.lons(1);
        
        
        % Check if reverse station pair has already been plotted
        stapairinv = [sta2,'_',sta1];
        if exist('existpair','var')
            if find(strncmp(stapairinv,existpair,length(stapairinv)))
                continue
            end
        end
        
        % Update some other useful variables
        dumsta2{nstapair} = sta2;
        npairall = npairall + 1; % number of total station pairs
        ccf_all{npairall} = ccf_filt{nstapair} ; % cell containing all ccf
        sta1sta2_dist_all(npairall) = sta1sta2_dist(nstapair); % vector containing distance between each station pair
        existpair(npairall) = {[sta1,'_',sta2]};
        
        % SNR
        [snr(npairall), signal_ind] = calc_SNR(ccf_filtered,min_grv,max_grv,sta1sta2_dist(nstapair),isfigure_snr);
        [snr_compare(npairall), ~] = calc_SNR(ccf_filtered_SNR,min_grv,max_grv,sta1sta2_dist(nstapair),isfigure_snr);
        dep1(npairall) = DEPTHS(strcmp(sta1,STAS));
        dep2(npairall) = DEPTHS(strcmp(sta2,STAS));
        
    end % ista2
    
    
    
    if IsFigure
        %----------- PLOT CCFs IN DISTANCE-TIME -------------%
        f101 = figure(101); clf; hold on;
        N= length(ccf_ifft);
        time = [-N/2:N/2]; % build lagtime vector for plotting
%         amp = 1e1; % amplitude to plot data
        indtime = find(abs(time)<=500); % Time index -500 to 500 seconds
        set(gca,'YDir','reverse');
        for istapair = 1: nstapair % loop over station pairs
            ccf_waveform = ccf_filt{istapair}(indtime(1):indtime(end)); % ccf at -500 to 500 seconds
            plot(time(indtime(1):indtime(end)),ccf_waveform*amp+sta1sta2_dist(istapair),'-k'); hold on;
%             text(0,stapairdist(istapair),dumsta2{istapair})
        end
        xlim([-500 500])
        xlabel('lag time (s)','fontsize',18,'fontweight','bold');
        ylabel('Distance (km)','fontsize',18,'fontweight','bold');
        title(['reference station:',sta1,'  filtered at ',num2str(coperiod(1)), ' - ',num2str(coperiod(2)),'(s)'],'fontsize',18,'fontweight','bold');
        set(gca,'fontsize',15);
        pause;
%         print(f101,'-dpdf',[figpath,'ccf',comp,'_',sta1,'.pdf']); % Save figure
    end
end % ista1


%% %----------- PLOT ALL CCFs STATION PAIRS IN DISTANCE-TIME -------------%
N= length(ccf_ifft);
time = [0:N-1]-floor(N/2);
time = [time(time<0), time(time>=0)];
% amp = 1e1;
if isploth20 && comp(1) == 'Z'
    amp = amp*1.5;
end
indtime_pos = find(abs(time)<=xlims(2) & time>=0);
indtime_neg = find(abs(time)<=xlims(2) & time<=0);

f102 = figure(102);
set(gcf, 'Color', 'w');
clf
hold on;
set(gca,'YDir','reverse');
dists = 0;
itrace = 0;
for istapair = 1: npairall
    if min(abs(sta1sta2_dist_all(istapair)-dists)) > trace_space && snr_compare(istapair) > snr_thresh ...
            && dep1(istapair) <= dep_tol(1) && dep2(istapair) <= dep_tol(2)
        itrace = itrace+1;
        dists(itrace) = sta1sta2_dist_all(istapair);
        
        % Plot positive
        ccf_waveform_all = ccf_all{istapair}(indtime_pos(1):indtime_pos(end)) / max(abs(ccf_all{istapair}(indtime_pos(1):indtime_pos(end))));
        plot(time(indtime_pos(1):indtime_pos(end)),ccf_waveform_all*amp+sta1sta2_dist_all(istapair),'-k','linewidth',1); hold on;
        
        % Plot negative
        ccf_waveform_all = ccf_all{istapair}(indtime_neg(1):indtime_neg(end)) / max(abs(ccf_all{istapair}(indtime_neg(1):indtime_neg(end))));
        plot(abs(time(indtime_neg(1):indtime_neg(end))),ccf_waveform_all*amp+sta1sta2_dist_all(istapair),'-r','linewidth',1); hold on;
    
    end
end
% xlim([-500 500])
xlim(xlims)
ylim(ylims);
xlabel('lag time (s)','fontsize',18);
ylabel('Distance (km)','fontsize',18);
title([comp(1),' : All non-repeated pairs, filtered at ',num2str(coperiod(1)), ' -',num2str(coperiod(2)),'(s)'],'fontsize',18);
set(gca,'fontsize',15);

% Plot Velocities
if isplotwin
    % Branches
    plot([min(sta1sta2_dist_all) max(sta1sta2_dist_all)]/max_grv,[min(sta1sta2_dist_all) max(sta1sta2_dist_all)],'color',[1 0 0],'linewidth',2);
    plot([min(sta1sta2_dist_all) max(sta1sta2_dist_all)]/-max_grv,[min(sta1sta2_dist_all) max(sta1sta2_dist_all)],'color',[1 0 0],'linewidth',2);
    plot([min(sta1sta2_dist_all) max(sta1sta2_dist_all)]/min_grv,[min(sta1sta2_dist_all) max(sta1sta2_dist_all)],'color',[1 0 0],'linewidth',2);
    plot([min(sta1sta2_dist_all) max(sta1sta2_dist_all)]/-min_grv,[min(sta1sta2_dist_all) max(sta1sta2_dist_all)],'color',[1 0 0],'linewidth',2);    
end

if isploth20 && comp(1) == 'Z'
    plot([min(sta1sta2_dist_all) max(sta1sta2_dist_all)]/h20_grv,[min(sta1sta2_dist_all) max(sta1sta2_dist_all)],'--','color',[0.5 0.5 1],'linewidth',2);
    plot([min(sta1sta2_dist_all) max(sta1sta2_dist_all)]/-h20_grv,[min(sta1sta2_dist_all) max(sta1sta2_dist_all)],'--','color',[0.5 0.5 1],'linewidth',2);
end

%pause;
%% Plot SNR Values
figure(101); clf;
plot([0 1000],[1 1],'-k','linewidth',3); hold on;
% plot(sta1sta2_dist_all,snr,'ok','linewidth',1,'MarkerFaceColor',[0.5 0.5 0.5],'markersize',8); hold on;
scatter(sta1sta2_dist_all,snr,80,mean([dep1' dep2'],2),'filled','MarkerEdgeColor',[0 0 0],'linewidth',1);
set(gca,'fontsize',16,'linewidth',2,'YScale','log');
grid on;
cb = colorbar;
ylabel(cb,'Average Depth (km)','fontsize',16);
xlabel('Distance','FontWeight','bold');
ylabel('SNR','FontWeight','bold');
xlim([0 500]);
ylim([1e-1 1e3]);


save2pdf([figpath,'SNR_',comp,'_tukeyfilt_',num2str(min_grv),'_',num2str(max_grv),'.pdf'],101,1000);
%%
% print(f102,'-dpdf',[figpath,'all_ccf',comp,'_tukeyfilt.pdf']); % Save figure
if isplotwin
    figname = [figpath,'all_ccf',comp,'_tukeyfilt_win',num2str(min_grv),'_',num2str(max_grv),'_TEI19.pdf'];
else
    figname = [figpath,'all_ccf',comp,'_tukeyfilt_TEI19.pdf'];
end
save2pdf(figname,f102,1000);
% export_fig(figname,'-pdf','-q100','-p0.02','-painters',f102)
% print(f102,'-dpdf',figname);
