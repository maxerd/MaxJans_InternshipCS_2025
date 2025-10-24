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

%%
clear MC

nMC = 20;

for i=1:nMC
    disp(['Iteration ',num2str(i)])
    try
        seperate_OLsim_sysID_V4
    catch
        disp('Problem in current iteration, trying again')
        seperate_OLsim_sysID_V4
    end

    MC.firstApprox{i} = sysID.par.firstAprrox.OL.raw;
    MC.initSys{i}     = sysID.par.initSys.sim.OL.raw;
    MC.fixedOrder{i}  = sysID.par.fixedOrder.sim.OL.raw;
    MC.straight{i}    = sysID.par.straight.sim.OL.raw;
end

%%
figure(201);clf
for i=1:nMC
    bode(MC.firstApprox{i});hold on
end
title('First order approximation, monte carlo');
grid minor

figure(202);clf
for i=1:nMC
    bode(MC.initSys{i});hold on
end
title('ssest using initial system, monte carlo');
grid minor

figure(203);clf
for i=1:nMC
    bode(MC.fixedOrder{i});hold on
end
title('ssest using a fixed order, monte carlo');
grid minor

figure(204);clf
for i=1:nMC
    bode(MC.straight{i});hold on
end
title('Straight data, monte carlo');
grid minor





























