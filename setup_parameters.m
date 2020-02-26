% SETUP PARAMETERS FOR RADON TRANSFORM
addpath('./functions/'); addpath('./functions/CG_methods/');

% Input data
comp = 'PP'; %'ZZ'; 'RR'; 'TT';
ccf = 'ccf_raw';
windir = 'window3hr';
% windir = 'window3hr_Zcorr_tiltcomp'; % Tilt & compliance corrected

% comp = 'ZZ';
% ccf = 'ccf_OBNW';
% windir = 'window3hr';

%% PARAMETERS
% Normalization option for plotting
is_globnorm = 0; % 1 for normalize radon panel by global max; 0 for column norm

% Define variables for calculating the LRT
maxiter = 10; % Maximum number of iterations
rthresh = 1e-6;
method = 'CGG_weight';
% method = 'CG_IRLS';
f_min = 1/40;%1/150;
f_max = 1/3;%1/20;
v_min = 1; %4;
v_max = 6; %8;
P_axis = [111/(v_max*1.1) : 0.1 : 111/(v_min*0.9)]; % s/deg
P_axis = P_axis / 111; %(s/km);

% Parameters for tracing dispersion curves
min_peak_prom = 0.5; % Minimum peak prominence, threshold for peak height
min_peak_dist = 1; % Minimum separation between chosen peaks [km/s]
Npers = 40; % Number for periods
% pers = logspace(log10(3),log10(40),Npers); % period vector
pers = 1./flip(linspace(f_min, f_max ,Npers));

%% SETUP PATHS
% Path to seismogram data structure
ndata = ['./data/',ccf,'/',windir,'/noise_',comp,'.mat'];

% Output LRT panel path (a1)
LRTmatpath = ['./LRT_mats/',ccf,'/',windir,'/',num2str(1/f_max),'_',num2str(1/f_min),'s/'];

% Output path for dispersion picks (a3)
picks_out_path = ['./LRT_picks/',ccf,'/',windir,'/',num2str(1/f_max),'_',num2str(1/f_min),'s/'];

% Figure path
figpath = ['./figs/',ccf,'/',windir,'/'];

%% LOAD DISPERSION
if strcmp(comp(1),'T')
    qfile = './qfiles/Nomelt_taper_eta_crust_INVpconstr_xi1.06_GRL19_ORCAiso_INV.t0to200.q';
else
    qfile = './qfiles/Nomelt_taper_eta_crust_INVpconstr_xi1.06_GRL19_ORCAiso_INV.s0to200.q';
end

BRANCHES=2;
for ii = 1:BRANCHES
    dat{ii} = readMINEOS_qfile_allper(qfile,ii-1);
    DISP(ii).n = ii-1;
    DISP(ii).cv =  dat{ii}.phv;
    DISP(ii).gv =  dat{ii}.grv;
    DISP(ii).cvq = dat{ii}.phvq;
    DISP(ii).Tq =  dat{ii}.Tq;
    DISP(ii).T =   dat{ii}.T;
%     plot(Tq(1:10:end),cvq(1:10:end),'--','color',[.5 .5 .5],'linewidth',1);   
end