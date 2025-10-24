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

% Ambient (starting) temperature
T_amb = 23;

%% Make the system used for data generation
disp('Constructing system model')
main_lumpedSystem_water_V3

G = sysTMC(32,[3 6]);

%% Define the identification signal
maxAmplitude = 120; % [W] Maximum exitation signal
positiveOnly = 1;  % [-] Define whether the exitation signal can be negative

disp('Generating multisine for identification')
dist     = genMultisine(fs, tMeas*3600*fs, 1, maxAmplitude,     positiveOnly)+25;

for i=1:20

seperate_CLsim_sysID_V5

G_CL_dir_MC(i) = G_CL_dir;
G_CL_cp_MC(i)  = G_CL_cp;
G_CL_2s_MC(i)  = G_CL_2s;

end

%%
figure(101);clf
for i=1:20
    bode(G_CL_dir_MC(i));hold on
    title('Direct method, monte carlo')
end
grid minor
figure(102);clf
for i=1:20
    bode(G_CL_cp_MC(i));hold on
    title('Indirect method, co-prime, monte carlo')
end
grid minor
figure(103);clf
for i=1:20
    bode(G_CL_2s_MC(i));hold on
    title('Indirect method, two stage, monte carlo')
end
grid minor





