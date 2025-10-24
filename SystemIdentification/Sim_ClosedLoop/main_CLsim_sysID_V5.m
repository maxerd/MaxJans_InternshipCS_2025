clear all
% close all
clc

%% Pre-requisites
CL_sim_setup

%% Main variables

% Figure number of of which all others are based
baseFig_sim = 8000;
plotAll = false;

% Save the identification?
saveID = true; % [true/false]

% Total simulation time
tMeas = 20; % [h]
% Sampling time of simulation
Ts = 10; % [s]
fs = 1/Ts; % [Hz]


%% Make the system used for data generation
disp('Constructing system model')
main_lumpedSystem_water_V3

G = sysTMC(32,[3 6]);

%% Define the identification signal
maxAmplitude = 120; % [W] Maximum exitation signal
positiveOnly = 1;  % [-] Define whether the exitation signal can be negative

disp('Generating multisine for identification')
dist     = genMultisine(fs, tMeas*3600*fs, 1, maxAmplitude,     positiveOnly)+25;

%% Run the closed loop simuation and identification
seperate_CLsim_sysID_V5

% Identification outputs
G_CL_dir;
G_CL_cp;
G_CL_2s;

%% Visualize the identfication and their respective validation
CL_sim_validate_all_param

%% CHECK SLIDES
% lecture 11, slide 26









