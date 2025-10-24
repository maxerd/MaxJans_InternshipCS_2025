
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
tMeas = 8; % [hours] Measurement time
% tMeas = 45; % [hours] Measurement time

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
maxAmplitude = 2; % [W] Maximum exitation signal
positiveOnly = 1;  % [-] Define whether the exitation signal can be negative

disp('Generating multisine for identification')
dist     = genMultisine(fs, tMeas*3600*fs, 1, maxAmplitude,     positiveOnly)+50;
dist_val = dist;

%% Run the open loop simuation and identification
compTrans = true;
newData   = true;
seperate_OLsim_sysID_V4

% sysID.nonPar.OL.frf.trd.raw;
% sysID.nonPar.OL.frf.lpm.raw;
% sysID.par.OL.time.firstApprox.raw;
% sysID.par.OL.time.initSys.sim.raw;
% sysID.par.OL.time.initSys.prd.raw;
% sysID.par.OL.time.fixedOrder.sim.raw;
% sysID.par.OL.time.fixedOrder.prd.raw;
% sysID.par.OL.freq.trd.sim.raw;
% sysID.par.OL.freq.lpm.prd.raw;
% 
% sysID.nonPar.OL.frf.trd.filt;
% sysID.nonPar.OL.frf.lpm.filt;
% sysID.par.OL.time.firstApprox.filt;
% sysID.par.OL.time.firstApprox.filt;
% sysID.par.OL.time.initSys.sim.filt;
% sysID.par.OL.time.initSys.prd.filt;
% sysID.par.OL.time.fixedOrder.sim.filt;
% sysID.par.OL.time.fixedOrder.prd.filt;
% sysID.par.OL.freq.trd.sim.filt;
% sysID.par.OL.freq.lpm.prd.filt;
% 
% sysID.par.OL.time.firstApproxExp;

%% Run the open loop identification on measurement
compTrans = true;
newData   = true;
dataDir_OL = 'C:\Users\maxja\Documents\(4)School\Master\Q9_Internship\matlabFiles\measurements\processedData\p__CB_none__ST_ms__SM_24w__SA_12w__DT_250829__MD_8h__WT_no__DS_1000.mat'; % Open loop data used for identification
seperate_OLmeas_sysID_V6

% sysID.nonPar.OL.frf.trd.raw;
% sysID.nonPar.OL.frf.lpm.raw;
% sysID.par.OL.time.firstApprox.raw;
% sysID.par.OL.time.initSys.sim.raw;
% sysID.par.OL.time.initSys.prd.raw;
% sysID.par.OL.time.fixedOrder.sim.raw;
% sysID.par.OL.time.fixedOrder.prd.raw;
% sysID.par.OL.freq.trd.sim.raw;
% sysID.par.OL.freq.lpm.prd.raw;
% 
% sysID.nonPar.OL.frf.trd.filt;
% sysID.nonPar.OL.frf.lpm.filt;
% sysID.par.OL.time.firstApprox.filt;
% sysID.par.OL.time.firstApprox.filt;
% sysID.par.OL.time.initSys.sim.filt;
% sysID.par.OL.time.initSys.prd.filt;
% sysID.par.OL.time.fixedOrder.sim.filt;
% sysID.par.OL.time.fixedOrder.prd.filt;
% sysID.par.OL.freq.trd.sim.filt;
% sysID.par.OL.freq.lpm.prd.filt;
% 
% sysID.par.OL.time.firstApproxExp;

%% Run the Neural network identification on open loop data
saveDataToNN(sysID.OLdataFilt.train.y,sysID.OLdataFilt.train.u)

saveFileNN = '_IDdata\NN\lstm_model_OLdata.keras';
epochs = 100;
rmseTrain = trainNN(saveFileNN,epochs);

[predNN_OL,errNN_OL] = predNN(saveFileNN,sysID.OLdataFilt.valid.y,sysID.OLdataFilt.valid.u);

figure(901);clf
    subplot(211)
        plot(predNN_OL);hold on
        plot(sysID.OLdataFilt.valid.y);hold on
    subplot(212)
        plot(errNN_OL);hold on
            legend([num2str(rms(errNN_OL))])

%% Run the Neural network training to see convergence time to measurement time
clear rmseTrain

tMeas_NN = 2:2:40; %[h]
% tMeas_NN = 2:0.5:7.5; %[h]
tMeas_NN_s = tMeas_NN*3600; %[s]
tMeas_NN_sample = tMeas_NN_s*fs; %[-]

saveFileNN = '_IDdata\NN\lstm_model_OLdata.keras';
epochs = 100;

for i=1:12
    saveDataToNN(sysID.OLdataFilt.train.y(1:tMeas_NN_sample(i)),sysID.OLdataFilt.train.u(1:tMeas_NN_sample(i)))
    
    rmseTrain{i} = trainNN(saveFileNN,epochs);
end

%%
n_trains = 12;
clear convEpoch


for i=1:n_trains
    rmseTrainSurf(i,:) = [rmseTrain{i} min(rmseTrain{i}).*ones(1,epochs-size(rmseTrain{i},2))];
    convEpoch(i) = size(rmseTrain{i},2);
end
trainData.val_rmse = rmseTrainSurf;
trainData.epochs = 1:epochs;
trainData.measTime = tMeas_NN(1:n_trains);

bestTrain = find(trainData.val_rmse(1:n_trains,end)==min(trainData.val_rmse(1:n_trains,end)),1,'first');

figure(10);
set(gcf,'Position',[600 300 900 400])
sgtitle('Analysis of prediction RMSE of Neural Network on validation data')
    subplot(3,2,[1 3 5])
        surf(trainData.measTime,(trainData.epochs),((trainData.val_rmse(1:n_trains,:)))')
        % set(gca,'YScale','log');
        set(gca,'ZScale','log');
        set(gca,'ColorScale','log');
        view([1 1 0.7])
            xlabel('Total time of training data [hours]')
            ylabel('Training epoch')
            zlabel('RMSE [degC]')
            title('Prediction RMSE over data length and training epoch')
    subplot(3,2,2)
        semilogy(trainData.epochs,((trainData.val_rmse(bestTrain,:))),LineWidth=lw);grid minor
            xlabel('Training epoch')
            ylabel('RMSE [degC]')
            title('Prediction RMSE of training with best final RMSE')
    ax324 = subplot(3,2,4);
        semilogy(trainData.measTime,((trainData.val_rmse(1:n_trains,end))),LineWidth=lw);grid minor
        % plot(trainData.measTime,((trainData.val_rmse(1:n_trains,end))));grid minor
            xlabel('Total time of training data available [hours]')
            ylabel('Best RMSE [degC]')
            title('Best RMSE for available time')
    ax326 = subplot(3,2,6);
        plot(trainData.measTime,convEpoch,LineWidth=lw);grid minor
            xlabel('Total time of training data available [hours]')
            ylabel('Epoch')
            title('Epochs needed to converge to a constant RMSE')
    linkaxes([ax324 ax326],'x')
    print('./_Figures/NN/ConvergencePlots_measurement','-depsc')

%% Make all the open loop figures for report/validation
valData = sysID.OLdataFilt.valid;

preds_meas = compare(valData.id,...
                    sysID.par.OL.time.firstApprox.filt,...
                    sysID.par.OL.time.fixedOrder1.prd.filt,...
                    sysID.par.OL.time.fixedOrder5.prd.filt,...
                    sysID.par.OL.time.initSys.prd.filt,...
                    sysID.par.OL.time.straight.prd.filt,1);
[predNN_OL,errNN_OL] = predNN(saveFileNN,valData.y,valData.u);

Ndelay = 1;

figure(17);clf;hold on;grid minor
set(gcf,'Position',[300 100 800 450])
    subplot(211);hold on;grid minor
        plot(sysID.OLdata.valid.tVec,predNN_OL)
        plot(sysID.OLdata.valid.tVec,(preds_meas{1}.y))
        plot(sysID.OLdata.valid.tVec,(preds_meas{2}.y))
        plot(sysID.OLdata.valid.tVec,(preds_meas{3}.y))
        plot(sysID.OLdata.valid.tVec,(preds_meas{4}.y))
        plot(sysID.OLdata.valid.tVec,(preds_meas{5}.y))
        plot(sysID.OLdata.valid.tVec,valData.y,'g--')
            xlabel('Time [s]')
            ylabel('Temperature [degC]')
            title('1-step ahead prediction of traditional and NN model')
            legend('Neural Network','First order approx.','Fixed order, nx=1','Fixed order, nx=5','Initialized system','Straight','Measurement data')
            % legend('Neural Network','Fixed order, nx=1','Fixed order, nx=5','Initialized system','Straight','Measurement data')

    subplot(212);hold on;grid minor
        plot(sysID.OLdata.valid.tVec,errNN_OL)
        plot(sysID.OLdata.valid.tVec,valData.y-(preds_meas{1}.y))
        plot(sysID.OLdata.valid.tVec,valData.y-(preds_meas{2}.y))
        plot(sysID.OLdata.valid.tVec,valData.y-(preds_meas{3}.y))
        plot(sysID.OLdata.valid.tVec,valData.y-(preds_meas{4}.y))
        plot(sysID.OLdata.valid.tVec,valData.y-(preds_meas{5}.y))
            xlabel('Time [s]')
            ylabel('Temperature error [degC]')
            title('1-step ahead prediction errors of traditional and NN model')
            legend(['Neural Network, RMSE: ',num2str(rms(errNN_OL))],...
                   ['First order approx, RMSE: ',num2str(rms(valData.y-(preds_meas{1}.y)))],...
                   ['Fixed order, nx=1, RMSE: ', num2str(rms(valData.y-(preds_meas{2}.y)))],...
                   ['Fixed order, nx=5, RMSE: ', num2str(rms(valData.y-(preds_meas{3}.y)))],...
                   ['Initialized system, RMSE: ',num2str(rms(valData.y-(preds_meas{4}.y)))],...
                   ['Straight, RMSE: ',          num2str(rms(valData.y-(preds_meas{5}.y)))])
            % legend(['Neural Network, RMSE: ',num2str(rms(errNN_OL))],...
            %        ['Fixed order, nx=1, RMSE: ' ,num2str(rms(valData.y-(preds_meas{2}.y)))],...
            %        ['Fixed order, nx=5, RMSE: ' ,num2str(rms(valData.y-(preds_meas{3}.y)))],...
            %        ['Initialized system, RMSE: ',num2str(rms(valData.y-(preds_meas{4}.y)))],...
            %        ['Straight, RMSE: ',          num2str(rms(valData.y-(preds_meas{5}.y)))])

    %%
firstOrderApprox
fixedOrderComparisons
timeApproachesComparisons
freqApproachesComparisons

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

%% Run the closed loop simuation and identification
% for i=1:10
seperate_CLmeas_sysID_V6

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

%% Run the Neural network identification on closed loop data
saveDataToNN(sysID.CLdataFilt.train.y,sysID.CLdataFilt.train.u)

saveFileNN = '_IDdata\NN\lstm_model_CLdata.keras';
epochs = 100;
trainNN(saveFileNN,epochs)

[predNN_CL,errNN_CL] = predNN(saveFileNN,sysID.CLdataFilt.valid.y,sysID.CLdataFilt.valid.u);
% predNN_CL = predNN(saveFileNN,sysID.dataFilt.OL.valid.y,sysID.dataFilt.OL.valid.u);

figure(902);clf
    subplot(211)
        plot(predNN_CL);hold on
        plot(sysID.CLdataFilt.valid.y);hold on
    subplot(212)
        plot(errNN_CL);hold on
            legend([num2str(rms(errNN_CL))])

%% Run the Neural network training to see convergence time to measurement time
clear rmseTrain

tMeas_NN = 2:2:34; %[h]
% tMeas_NN = 2.5:1:6.5; %[h]
tMeas_NN_s = tMeas_NN*3600; %[s]
tMeas_NN_sample = tMeas_NN_s*fs; %[-]

saveFileNN = '_IDdata\NN\lstm_model_CLdata.keras';
epochs = 100;
n_trains = size(tMeas_NN,2);

for i=1:n_trains
    saveDataToNN(sysID.CLdataFilt.train.y(1:tMeas_NN_sample(i)),sysID.CLdataFilt.train.u(1:tMeas_NN_sample(i)))
    
    rmseTrain{i} = trainNN(saveFileNN,epochs);
end

%%
clear convEpoch rmseTrainSurf


for i=1:n_trains
    rmseTrainSurf(i,:) = [rmseTrain{i} min(rmseTrain{i}).*ones(1,epochs-size(rmseTrain{i},2))];
    convEpoch(i) = size(rmseTrain{i},2);
end
trainData.val_rmse = rmseTrainSurf;
trainData.epochs = 1:epochs;
trainData.measTime = tMeas_NN(1:n_trains);

bestTrain = find(trainData.val_rmse(1:n_trains,end)==min(trainData.val_rmse(1:n_trains,end)),1,'first');

figure(10);
set(gcf,'Position',[600 300 900 400])
sgtitle('Analysis of prediction RMSE of Neural Network on validation data')
    subplot(3,2,[1 3 5])
        surf(trainData.measTime,(trainData.epochs),((trainData.val_rmse(1:n_trains,:)))')
        % set(gca,'YScale','log');
        set(gca,'ZScale','log');
        set(gca,'ColorScale','log');
        view([1 1 0.7])
            xlabel('Total time of training data [hours]')
            ylabel('Training epoch')
            zlabel('RMSE [degC]')
            title('Prediction RMSE over data length and training epoch')
    subplot(3,2,2)
        semilogy(trainData.epochs,((trainData.val_rmse(bestTrain,:))),LineWidth=lw);grid minor
            xlabel('Training epoch')
            ylabel('RMSE [degC]')
            title('Prediction RMSE of training with best final RMSE')
    ax324 = subplot(3,2,4);
        semilogy(trainData.measTime,((trainData.val_rmse(1:n_trains,end))),LineWidth=lw);grid minor
        % plot(trainData.measTime,((trainData.val_rmse(1:n_trains,end))));grid minor
            xlabel('Total time of training data available [hours]')
            ylabel('Best RMSE [degC]')
            title('Best RMSE for available time')
    ax326 = subplot(3,2,6);
        plot(trainData.measTime,convEpoch,LineWidth=lw);grid minor
            xlabel('Total time of training data available [hours]')
            ylabel('Epoch')
            title('Epochs needed to converge to a constant RMSE')
    linkaxes([ax324 ax326],'x')
    print('./_Figures/NN/ConvergencePlots_measurement','-depsc')




%% Compare the results from OL and CL identification

figure(51);clf
    bode(sysID.par.OL.time.firstApprox.filt,sysID.par.OL.time.initSys.prd.filt,sysID.par.OL.time.fixedOrder.prd.filt,sysID.par.OL.time.straight.prd.filt,...
         sysID.par.CL.time.direct.filt,sysID.par.CL.time.coprime.filt,sysID.par.CL.time.twostage.filt,G(1,1),'g--');grid minor
    xlim([1e-5 1])
        title('Comparison of identification methods')
        legend('OL, first order approximation','OL, initialized system','OL, fixed order','OL, straight data','CL, direct method','CL, co-prime method','CL, two-stage method','location','best')
figure(52);clf
    bode(sysID.par.OL.time.firstApprox.filt,sysID.par.OL.time.initSys.prd.filt,sysID.par.OL.time.fixedOrder.prd.filt,sysID.par.OL.time.straight.prd.filt,...
         sysID.par.CL.time.direct.filt,sysID.par.CL.time.coprime.filt,sysID.par.CL.time.twostage.filt,G(1,1),'g--');grid minor
    xlim([1e-5 1])
        title('Comparison of identification methods')
        legend('OL, first order approximation','OL, initialized system','OL, fixed order','OL, straight data','CL, direct method','CL, two-stage method','location','best')


figure(53);clf
    simpleBodemag(sysID.par.OL.time.firstApprox.filt,'Hz',bodeRange);hold on;grid minor
    simpleBodemag(sysID.par.OL.time.initSys.prd.filt,'Hz',bodeRange)
    simpleBodemag(sysID.par.OL.time.fixedOrder.prd.filt,'Hz',bodeRange)
    simpleBodemag(sysID.par.OL.time.straight.prd.filt,'Hz',bodeRange)
    simpleBodemag(sysID.par.CL.time.direct.filt,'Hz',bodeRange)
    simpleBodemag(sysID.par.CL.time.coprime.filt,'Hz',bodeRange)
    simpleBodemag(sysID.par.CL.time.twostage.filt,'Hz',bodeRange)
    simpleBodemag(G(1,1),'Hz',bodeRange,'g--')
    xlim([1e-5 1])
        title('Comparison of identification methods')
        legend('OL, first order approximation','OL, initialized system','OL, fixed order','OL, straight data','CL, direct method','CL, co-prime method','CL, two-stage method','location','best')

%%
figure(61);clf
    compare(sysID.OLdataFilt.valid.id,sysID.par.OL.time.firstApprox.filt,sysID.par.OL.time.initSys.prd.filt,sysID.par.OL.time.fixedOrder.prd.filt,sysID.par.OL.time.straight.prd.filt,...
         d2c(sysID.par.CL.time.direct.filt),d2c(sysID.par.CL.time.coprime.filt),d2c(sysID.par.CL.time.twostage.filt),G(1,1),'g--',inf);grid minor
        title('Comparison of identification methods of the TMC setup')
        legend('Data','OL, first order approximation','OL, initialized system','OL, fixed order','OL, straight data','CL, direct method','CL, co-prime method','CL, two-stage method','Model','location','best')
figure(62);clf
    compare(sysID.CLdataFilt.valid.direct,sysID.par.OL.time.firstApprox.filt,sysID.par.OL.time.initSys.prd.filt,sysID.par.OL.time.fixedOrder.prd.filt,sysID.par.OL.time.straight.prd.filt,...
         sysID.par.CL.time.direct.filt,sysID.par.CL.time.coprime.filt,sysID.par.CL.time.twostage.filt,G(1,1),'g--',inf);grid minor
        title('Simulation comparison of identification methods')
        legend('Data','OL, first order approximation','OL, initialized system','OL, fixed order','OL, straight data','CL, direct method','CL, co-prime method','CL, two-stage method','Model','location','best')
% figure(161);clf
%     compare(sysID.OLdataFilt.valid.id,sysID.par.OL.time.firstApprox.filt,sysID.par.OL.time.initSys.prd.filt,sysID.par.OL.time.fixedOrder.prd.filt,sysID.par.OL.time.straight.prd.filt,...
%          sysID.par.CL.time.direct.filt,sysID.par.CL.time.coprime.filt,sysID.par.CL.time.twostage.filt,G(1,1),'g--');grid minor
%         title('Comparison of identification methods of the TMC setup')
%         legend('Data','OL, first order approximation','OL, initialized system','OL, fixed order','OL, straight data','CL, direct method','CL, two-stage method','location','best')
% figure(162);clf
%     compare(sysID.CLdataFilt.valid.direct,sysID.par.OL.time.firstApprox.filt,sysID.par.OL.time.initSys.prd.filt,sysID.par.OL.time.fixedOrder.prd.filt,sysID.par.OL.time.straight.prd.filt,...
%          sysID.par.CL.time.direct.filt,sysID.par.CL.time.coprime.filt,sysID.par.CL.time.twostage.filt,G(1,1),'g--');grid minor
%         title('Simulation comparison of identification methods')
%         legend('Data','OL, first order approximation','OL, initialized system','OL, fixed order','OL, straight data','CL, direct method','CL, two-stage method','location','best')




















