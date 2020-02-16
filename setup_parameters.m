
%% SETUP PATHS
addpath('./functions/'); addpath('./functions/CG_methods/');
% Input data
ndata = './pa5_5km/Synth_120W_150W.mat';
% Output LRT panel path (a1)
LRTmatpath = './LRT_mats/';
% Output path for dispersion picks (a3)
out_path = './LRT_picks/';

%% PARAMETERS
% Normalization option for plotting
is_globnorm = 1; % 1 for normalize radon panel by global max; 0 for column norm

% Define variables for calculating the LRT
maxiter = 10; % Maximum number of iterations
rthresh = 1e-6;
method = 'CGG_weight';
% method = 'CG_IRLS';
f_min = 1/150;
f_max = 1/20;
v_min = 4;
v_max = 8;
P_axis = [111/(v_max*1.1) : 0.1 : 111/(v_min*0.9)]; % s/deg
P_axis = P_axis / 111; %(s/km);

% Parameters for tracing dispersion curves
min_peak_prom = 0.3; % Minimum peak prominence, threshold for peak height
min_peak_dist = 0.1; % Minimum separation between chosen peaks [km/s]
Npers = 25; % Number for periods
pers = logspace(log10(20),log10(150),Npers); % period vector 
