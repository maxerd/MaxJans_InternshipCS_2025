clear all
% close all
clc

%% Pre-requisites
OL_sim_setup

%% Main variables

% Figure number of of which all others are based
baseFig_sim = 8000;
plotAll = false;

% Save the identification?
saveID = false; % [true/false]

% Define some time variables
Ts = 10;  % [s]     Sampling time
fs = 1/Ts;     % [Hz]    Sampling frequency
tMeas = 8; % [hours] Measurement time

% Ambient (starting) temperature
T_amb = 23; % [degC] Ambient Temperature

%% Make the system used for data generation
disp('Constructing system model')
main_lumpedSystem_water_V3

G = sysTMC(32,[3 6]);

%% Define the input signal
maxAmplitude = 120; % [W] Maximum exitation signal
positiveOnly = 1;  % [-] Define whether the exitation signal can be negative

disp('Generating multisine for identification')
dist     = genMultisine(fs, tMeas*3600*fs, 1, maxAmplitude,     positiveOnly);
dist_val = dist;

%% Run the open loop simuation and identification
seperate_OLsim_sysID_V4

% Identification outputs
sysID.par.firstAprrox.OL.raw;
sysID.par.initSys.sim.OL.raw;
sysID.par.fixedOrder.sim.OL.raw;
sysID.par.straight.sim.OL.raw;

sysID.par.firstAprrox.OL.filt;
sysID.par.initSys.sim.OL.filt;
sysID.par.fixedOrder.sim.OL.filt;
sysID.par.straight.sim.OL.filt;















