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
comps = {'ZZ','PP'}; %{'ZZ','RR','PP'}; %'PZ'; %'PP'; %'ZZ'; %'RR'; %'TT';
amp = 8e0;
windir = 'window3hr';
% windir = 'window3hr_Zcorr_tiltcomp';


max_grv = inf; %5.5;
min_grv = 0.7; %1.6
xlims = [-500 500];
ylims = [0 450];

ccf_path = '../ccf/ccf_raw/';
% ccf_path = '/data/irma6/jrussel/YoungPacificORCA/ccf_FTN/';

data_dir = './data/ccf_raw/';

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

M = [];
Delta = [];
t = [];

for icomp = 1:length(comps)
    comp = comps{icomp};
    %%
    ccf_path = [ccf_stack_path,'ccf',comp,'/',];
    npairall = 0;
    %------------ LOAD DATA AND PLOT IN TIME DOMAIN -------------%
    for ista1=1:nsta % loop over all stations
        sta1=char(stalist(ista1,:));
        sta1dir=[ccf_path,sta1]; % dir to have all cross terms about this central station
        nstapair = 0;
        for ista2 = 1: nsta % loop over station pairs
            sta2 = char(stalist(ista2,:));

            % if same station, skip
            if(strcmp(sta1,sta2))
                continue
            end

            % Check if reverse station pair has already been plotted
            npairall = npairall+1;
            stapairinv = [sta2,'_',sta1];
            if exist('existpair','var')
                if find(strncmp(stapairinv,existpair,length(stapairinv)))
                    continue
                end
            end
            existpair(npairall) = {[sta1,'_',sta2]};

            filename = sprintf('%s/%s_%s_f.mat',sta1dir,sta1,sta2);

            if ~exist(filename,'file') % check that ccf file exists
                disp(['not exist ',filename])
                continue;
            end
            nstapair = nstapair + 1;

    %         % Want sta1 to be closest to the coast so waves at -lag travel
    %         % towards coast and waves at +lag travel away from coast.
    %         filename_sta1sta2 = filename;
    %         filename_sta2sta1 = [ccf_path,sta2,'/',sta2,'_',sta1,'_f.mat'];
    %         test = load(filename_sta1sta2);
    %         sta1lon = test.stapairsinfo.lons(1);
    %         sta2lon = test.stapairsinfo.lons(2);
    %         if sta1lon > sta2lon
    %             filename = filename_sta2sta1;
    %         else
    %             filename = filename_sta1sta2;
    %         end

            %----------- LOAD DATA -------------%
            data = load(filename);
            ccf = data.coh_sum./data.coh_num;
            ccf(isnan(ccf)) = 0;
            if size(ccf,1) == 1
                ccf = ccf';
            end

            %%
            %----------- Frequency ==> Time domain -------------%
            r = deg2km(distance(data.stapairsinfo.lats(1),data.stapairsinfo.lons(1),data.stapairsinfo.lats(2),data.stapairsinfo.lons(2)));
            N = length(ccf);
            ccf_ifft = real(ifft(ccf,N)); % inverse FFT to get time domain
            ccf_ifft = fftshift(ccf_ifft); % rearrange values as [-lag lag]
            ccf_ifft = detrend(ccf_ifft);
            ccf_ifft = cos_taper(ccf_ifft);

            % Build time axis [-lag:0:lag] and index positive and negative waveforms
            time = [0:N-1]-floor(N/2);
            time = [time(time<0), time(time>=0)];
            indtime_pos = find(abs(time)<=xlims(2) & time>=0);
            indtime_neg = find(abs(time)<=xlims(2) & time<=0);
            ccf_ifft_pos = cos_taper(detrend(ccf_ifft(indtime_pos)));
            ccf_ifft_neg = flip(cos_taper(detrend(ccf_ifft(indtime_neg))));

    %         % flip -lag to positive +lag
    %         M = [M; ccf_ifft_pos; ccf_ifft_neg ];
    %         Delta = [Delta; r; r ];
    %         t = time(indtime_pos);

            % average +lag and -lag
            M = [M; mean([ccf_ifft_pos; ccf_ifft_neg],1)];
            Delta = [Delta; r ];
            t = time(indtime_pos);

    %         % Keep -lag:+lag
    %         M = [M; cos_taper(detrend(ccf_ifft(abs(time)<=xlims(2)))) ];
    %         Delta = [Delta; r ];
    %         t = time(abs(time)<=xlims(2));


        end % ista2
    end % ista1
end

% Sort data by distance
[Delta, I_srt] = sort(Delta);
M = M(I_srt,:);

%% %----------- PLOT ALL CCFs STATION PAIRS IN DISTANCE-TIME -------------%
f102 = figure(102);
set(gcf, 'Color', 'w');
clf
hold on;
set(gca,'YDir','reverse');
dists = 0;
itrace = 0;
for istapair = 1: size(M,1)
    itrace = itrace+1;
    ccf_waveform_all = M(istapair,:) / max(abs(M(istapair,:)));
    plot(t,ccf_waveform_all*amp+Delta(istapair),'-k','linewidth',1); hold on;
end
xlim([min(t) max(t)])
% xlim([0 max(xlims)])
ylim(ylims);
xlabel('lag time (s)','fontsize',18);
ylabel('Distance (km)','fontsize',18);
set(gca,'fontsize',15);


% save2pdf([figpath,'SNR_',comp,'_tukeyfilt_',num2str(min_grv),'_',num2str(max_grv),'.pdf'],101,1000);

if ~exist(data_dir)
    mkdir(data_dir)
end
save([data_dir,'noise_',cell2mat(comps),'.mat'],'M','Delta','t');