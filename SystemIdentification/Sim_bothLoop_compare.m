clear all
% close all
clc

%% Pre-requisites
all_sim_setup

%% Main variables

% Figure number of of which all others are based
baseFig_sim = 9000;
plotAll = false;

% Save the identification?
saveID = false; % [true/false]

% Define some time variables
Ts = 10;  % [s]     Sampling time
fs = 1/Ts;     % [Hz]    Sampling frequency
tMeas = 20; % [hours] Measurement time

% Ambient (starting) temperature
T_amb = 23; % [degC] Ambient Temperature

%% Make the system used for data generation
disp('Constructing system model')
main_lumpedSystem_water_V3

G = sysTMC(32,[3 6]);

%% Define the input signal for open loop
maxAmplitude = 120; % [W] Maximum exitation signal
positiveOnly = 1;  % [-] Define whether the exitation signal can be negative

disp('Generating multisine for identification')
dist     = genMultisine(fs, tMeas*3600*fs, 1, maxAmplitude,     positiveOnly);
dist_val = dist;

%% Define the input signal for closed loop
% maxAmplitude = 3; % [W] Maximum exitation signal
% positiveOnly = 1;  % [-] Define whether the exitation signal can be negative
% 
% disp('Generating multisine for identification')
% dist     = genMultisine(fs, tMeas*3600*fs, 1, maxAmplitude,     positiveOnly)+80;
% dist_val = dist;

%% Run the open loop simuation and identification
compTrans = true;
newData   = false;
seperate_OLsim_sysID_V4

sysID.nonPar.OL.frf.trd.raw;
sysID.nonPar.OL.frf.lpm.raw;
sysID.par.OL.time.firstApprox.raw;
sysID.par.OL.time.initSys.sim.raw;
sysID.par.OL.time.initSys.prd.raw;
sysID.par.OL.time.fixedOrder.sim.raw;
sysID.par.OL.time.fixedOrder.prd.raw;
sysID.par.OL.freq.trd.sim.raw;
sysID.par.OL.freq.lpm.prd.raw;

sysID.nonPar.OL.frf.trd.filt;
sysID.nonPar.OL.frf.lpm.filt;
sysID.par.OL.time.firstApprox.filt;
sysID.par.OL.time.firstApprox.filt;
sysID.par.OL.time.initSys.sim.filt;
sysID.par.OL.time.initSys.prd.filt;
sysID.par.OL.time.fixedOrder.sim.filt;
sysID.par.OL.time.fixedOrder.prd.filt;
sysID.par.OL.freq.trd.sim.filt;
sysID.par.OL.freq.lpm.prd.filt;

sysID.par.OL.time.firstApproxExp;


%% Run the Neural network identification on open loop data
saveDataToNN(sysID.OLdataFilt.train.y,sysID.OLdataFilt.train.u)

saveFileNN = '_IDdata\NN\lstm_model_OLdata.keras';
trainNN(saveFileNN)

[predNN_OL,errNN_OL] = predNN(saveFileNN,sysID.OLdataFilt.valid.y,sysID.OLdataFilt.valid.u);

figure(901);clf
    subplot(211)
        plot(predNN_OL);hold on
        plot(sysID.OLdataFilt.valid.y);hold on
    subplot(212)
        plot(errNN_OL);hold on
            legend([num2str(rms(errNN_OL))])

%% Run the closed loop simuation and identification
% for i=1:10
seperate_CLsim_sysID_V5

% Identification outputs
sysID.nonPar.CL.frf.direct.raw;
sysID.nonPar.CL.frf.classic.raw;
sysID.nonPar.CL.frf.coprime.raw;
sysID.par.CL.time.direct.raw;
sysID.par.CL.time.classic.raw;
sysID.par.CL.time.coprime.raw;
sysID.par.CL.time.twostage.raw;

sysID.nonPar.CL.frf.direct.filt;
sysID.nonPar.CL.frf.classic.filt;
sysID.nonPar.CL.frf.coprime.filt;
sysID.par.CL.time.direct.filt;
sysID.par.CL.time.classic.filt;
sysID.par.CL.time.coprime.filt;
sysID.par.CL.time.twostage.filt;

% pred_dir = compare(sysID.data.CL.valid.direct,G_CL_dir);
% err_rms_dir(i) = rmse(pred_dir.y,sysID.data.CL.valid.y);
% pred_cp = compare(sysID.data.CL.valid.direct,G_CL_cp);
% err_rms_cp(i) = rmse(pred_cp.y,sysID.data.CL.valid.y);
% pred_2s = compare(sysID.data.CL.valid.direct,G_CL_2s);
% err_rms_2s(i) = rmse(pred_2s.y,sysID.data.CL.valid.y);
% end

%% Run the Neural network identification on closed loop data
saveDataToNN(sysID.CLdataFilt.train.y,sysID.CLdataFilt.train.u)

saveFileNN = '_IDdata\NN\lstm_model_CLdata.keras';
trainNN(saveFileNN)

[predNN_CL,errNN_CL] = predNN(saveFileNN,sysID.CLdataFilt.valid.y,sysID.CLdataFilt.valid.u);
% predNN_CL = predNN(saveFileNN,sysID.dataFilt.OL.valid.y,sysID.dataFilt.OL.valid.u);

figure(902);clf
    subplot(211)
        plot(predNN_CL);hold on
        plot(sysID.CLdataFilt.valid.y);hold on
    subplot(212)
        plot(errNN_CL);hold on
            legend([num2str(rms(errNN_CL))])

%% Make all the figures for report/validation
firstOrderApprox
fixedOrderComparisons
timeApproachesComparisons
freqApproachesComparisons




%% Compare the results from OL and CL identification

figure(51);clf
    bode(sysID.par.firstAprrox.OL.filt,sysID.par.initSys.sim.OL.filt,sysID.par.fixedOrder.sim.OL.filt,sysID.par.straight.sim.OL.filt,...
         G_CL_dir,G_CL_cp,G_CL_2s,G(1,1),'g--');grid minor
    xlim([1e-5 1])
        title('Comparison of identification methods')
        legend('OL, first order approximation','OL, initialized system','OL, fixed order','OL, straight data','CL, direct method','CL, co-prime method','CL, two-stage method','location','best')
figure(52);clf
    bode(sysID.par.firstAprrox.OL.filt,sysID.par.initSys.sim.OL.filt,sysID.par.fixedOrder.sim.OL.filt,sysID.par.straight.sim.OL.filt,...
         G_CL_dir,G_CL_2s,G(1,1),'g--');grid minor
    xlim([1e-5 1])
        title('Comparison of identification methods')
        legend('OL, first order approximation','OL, initialized system','OL, fixed order','OL, straight data','CL, direct method','CL, two-stage method','location','best')

figure(61);clf
    compare(sysID.OLdataFilt.valid.id,sysID.par.firstAprrox.OL.filt,sysID.par.initSys.sim.OL.filt,sysID.par.fixedOrder.sim.OL.filt,sysID.par.straight.sim.OL.filt,...
         G_CL_dir,G_CL_cp,G_CL_2s,G(1,1),'g--');grid minor
        title('Comparison of identification methods of the TMC setup')
        legend('Data','OL, first order approximation','OL, initialized system','OL, fixed order','OL, straight data','CL, direct method','CL, co-prime method','CL, two-stage method','location','best')
figure(62);clf
    compare(sysID.CLdataFilt.valid.direct,sysID.par.firstAprrox.OL.filt,sysID.par.initSys.sim.OL.filt,sysID.par.fixedOrder.sim.OL.filt,sysID.par.straight.sim.OL.filt,...
         G_CL_dir,G_CL_cp,G_CL_2s,G(1,1),'g--');grid minor
        title('Simulation comparison of identification methods')
        legend('Data','OL, first order approximation','OL, initialized system','OL, fixed order','OL, straight data','CL, direct method','CL, co-prime method','CL, two-stage method','location','best')
figure(161);clf
    compare(sysID.OLdataFilt.valid.id,sysID.par.firstAprrox.OL.filt,sysID.par.initSys.sim.OL.filt,sysID.par.fixedOrder.sim.OL.filt,sysID.par.straight.sim.OL.filt,...
         G_CL_dir,G_CL_2s,G(1,1),'g--');grid minor
        title('Comparison of identification methods of the TMC setup')
        legend('Data','OL, first order approximation','OL, initialized system','OL, fixed order','OL, straight data','CL, direct method','CL, two-stage method','location','best')
figure(162);clf
    compare(sysID.CLdataFilt.valid.direct,sysID.par.firstAprrox.OL.filt,sysID.par.initSys.sim.OL.filt,sysID.par.fixedOrder.sim.OL.filt,sysID.par.straight.sim.OL.filt,...
         G_CL_dir,G_CL_2s,G(1,1),'g--');grid minor
        title('Simulation comparison of identification methods')
        legend('Data','OL, first order approximation','OL, initialized system','OL, fixed order','OL, straight data','CL, direct method','CL, two-stage method','location','best')
























