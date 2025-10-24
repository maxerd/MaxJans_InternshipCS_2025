
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

N = tMeas*3600*fs; % [-] Measurement samples

% Ambient (starting) temperature
T_amb = 23; % [degC] Ambient Temperature

%% Make the system used for data generation
disp('Constructing system model')
main_lumpedSystem_water_V3

G = sysTMC(32,[3 6]);

v = 0.15.*randn(N,1);

%% Define the input signal for open loop
maxAmplitude = 120; % [W] Maximum exitation signal
positiveOnly = 1;  % [-] Define whether the exitation signal can be negative

disp('Generating multisine for identification')
dist     = genMultisine(fs, tMeas*3600*fs, 1, maxAmplitude,     positiveOnly);
dist_val = dist;

%% Run the open loop simuation and identification
compTrans = true;
newData   = true;
seperate_OLsim_sysID_V4


%% Define the input signal for closed loop
maxAmplitude = 20; % [W] Maximum exitation signal
positiveOnly = 1;  % [-] Define whether the exitation signal can be negative

disp('Generating multisine for identification')
dist     = genMultisine(fs, tMeas*3600*fs, 1, maxAmplitude,     positiveOnly)+118+25;
dist_val = dist;

%% Run the closed loop simuation and identification
compTrans = true;
newData   = true;
seperate_CLsim_sysID_V5_2

%% Compare FRF's from both open and closed loop
% The two different figures show the bodes from non-compensated and
% compensated data respectively

figure(3101);clf
set(gcf,'position',[300 100 1200 550])
    subplot(221)
        simpleBodemag(sysID.nonPar.OL.frf.trd.raw,'Hz','b',1);hold on;grid minor
        simpleBodemag(sysID.nonPar.CL.frf.trd.direct.raw,'Hz' ,'r',1)
        simpleBodemag(sysID.nonPar.CL.frf.trd.classic.raw,'Hz','k',1)
        simpleBodemag(sysID.nonPar.CL.frf.trd.coprime.raw,'Hz','m',1)
        simpleBodemag(G(1,1)-sysID.nonPar.OL.frf.trd.raw,'Hz','b--')
        simpleBodemag(G(1,1)-sysID.nonPar.CL.frf.trd.direct.raw,'Hz' ,'r--')
        simpleBodemag(G(1,1)-sysID.nonPar.CL.frf.trd.classic.raw,'Hz','k--')
        simpleBodemag(G(1,1)-sysID.nonPar.CL.frf.trd.coprime.raw,'Hz','m--')
        simpleBodemag(G(1,1),'Hz','g--')
        ylim([-100 25])
        xlim([0.0001 0.1])
            title('Bodes of non-parametric identifiations on simulation data, traditional method')
    subplot(223)
        simpleBodephase(sysID.nonPar.OL.frf.trd.raw,'Hz','wrap','b',1);hold on;grid minor
        simpleBodephase(sysID.nonPar.CL.frf.trd.direct.raw ,'Hz','wrap','r',1)
        simpleBodephase(sysID.nonPar.CL.frf.trd.classic.raw,'Hz','wrap','k',1)
        simpleBodephase(sysID.nonPar.CL.frf.trd.coprime.raw,'Hz','wrap','m',1)
        simpleBodephase(G(1,1),'Hz',logspace(-4,0,1000),'wrap','g--')
        xlim([0.0001 0.1])
            legend('OL, traditional','CL, Direct','CL, Indirect classic','CL, Indirect coprime','Data generating model')
    subplot(222)
        simpleBodemag(sysID.nonPar.OL.frf.lpm.raw,'Hz','b',1);hold on;grid minor
        simpleBodemag(sysID.nonPar.CL.frf.lpm.direct.raw ,'Hz','r',1)
        simpleBodemag(sysID.nonPar.CL.frf.lpm.classic.raw,'Hz','k',1)
        simpleBodemag(sysID.nonPar.CL.frf.lpm.coprime.raw,'Hz','m',1)
        simpleBodemag(G(1,1)-sysID.nonPar.OL.frfComp.lpm.raw,'Hz','b--')
        simpleBodemag(G(1,1)-sysID.nonPar.CL.frf.lpm.direct.raw ,'Hz','r--')
        simpleBodemag(G(1,1)-sysID.nonPar.CL.frf.lpm.classic.raw,'Hz','k--')
        simpleBodemag(G(1,1)-sysID.nonPar.CL.frf.lpm.coprime.raw,'Hz','m--')
        simpleBodemag(G(1,1),'Hz',logspace(-4,0,1000),'g--')
        ylim([-100 25])
        xlim([0.0001 0.1])
            title('Bodes of non-parametric identifiations on simulation data, LPM')
    subplot(224)
        simpleBodephase(sysID.nonPar.OL.frf.lpm.raw,'Hz','wrap','b',1);hold on;grid minor
        simpleBodephase(sysID.nonPar.CL.frf.lpm.direct.raw ,'Hz','wrap','r',1)
        simpleBodephase(sysID.nonPar.CL.frf.lpm.classic.raw,'Hz','wrap','k',1)
        simpleBodephase(sysID.nonPar.CL.frf.lpm.coprime.raw,'Hz','wrap','m',1)
        simpleBodephase(G(1,1),'Hz',logspace(-4,0,1000),'wrap','g--')
        xlim([0.0001 0.1])
            legend('OL, LPM','CL, Direct','CL, Indirect classic','CL, Indirect coprime','Data generating model')
    print('./_Figures/nonParam/nonCompensated_FRF_tradLPM_8h','-depsc')

figure(3102);clf
set(gcf,'position',[300 100 1200 550])
    subplot(221)
        simpleBodemag(sysID.nonPar.OL.frfComp.trd.raw,'Hz','b',1);hold on;grid minor
        simpleBodemag(sysID.nonPar.CL.frfComp.trd.direct.raw,'Hz' ,'r',1)
        simpleBodemag(sysID.nonPar.CL.frfComp.trd.classic.raw,'Hz','k',1)
        simpleBodemag(sysID.nonPar.CL.frfComp.trd.coprime.raw,'Hz','m',1)
        simpleBodemag(G(1,1)-sysID.nonPar.OL.frfComp.trd.raw,'Hz','b--')
        simpleBodemag(G(1,1)-sysID.nonPar.CL.frfComp.trd.direct.raw,'Hz' ,'r--')
        simpleBodemag(G(1,1)-sysID.nonPar.CL.frfComp.trd.classic.raw,'Hz','k--')
        simpleBodemag(G(1,1)-sysID.nonPar.CL.frfComp.trd.coprime.raw,'Hz','m--')
        simpleBodemag(G(1,1),'Hz','g--')
        ylim([-100 25])
        xlim([0.0001 0.1])
            title('Bodes of non-parametric identifiations on simulation data, traditional method')
    subplot(223)
        simpleBodephase(sysID.nonPar.OL.frfComp.trd.raw,'Hz','wrap','b',1);hold on;grid minor
        simpleBodephase(sysID.nonPar.CL.frfComp.trd.direct.raw ,'Hz','wrap','r',1)
        simpleBodephase(sysID.nonPar.CL.frfComp.trd.classic.raw,'Hz','wrap','k',1)
        simpleBodephase(sysID.nonPar.CL.frfComp.trd.coprime.raw,'Hz','wrap','m',1)
        simpleBodephase(G(1,1),'Hz',logspace(-4,0,1000),'wrap','g--')
        xlim([0.0001 0.1])
            legend('OL, traditional','CL, Direct','CL, Indirect classic','CL, Indirect coprime','Data generating model')
    subplot(222)
        simpleBodemag(sysID.nonPar.OL.frfComp.lpm.raw,'Hz','b',1);hold on;grid minor
        simpleBodemag(sysID.nonPar.CL.frfComp.lpm.direct.raw ,'Hz','r',1)
        simpleBodemag(sysID.nonPar.CL.frfComp.lpm.classic.raw,'Hz','k',1)
        simpleBodemag(sysID.nonPar.CL.frfComp.lpm.coprime.raw,'Hz','m',1)
        simpleBodemag(G(1,1)-sysID.nonPar.OL.frfComp.lpm.raw,'Hz','b--')
        simpleBodemag(G(1,1)-sysID.nonPar.CL.frfComp.lpm.direct.raw ,'Hz','r--')
        simpleBodemag(G(1,1)-sysID.nonPar.CL.frfComp.lpm.classic.raw,'Hz','k--')
        simpleBodemag(G(1,1)-sysID.nonPar.CL.frfComp.lpm.coprime.raw,'Hz','m--')
        simpleBodemag(G(1,1),'Hz',logspace(-4,0,1000),'g--')
        ylim([-100 25])
        xlim([0.0001 0.1])
            title('Bodes of non-parametric identifiations on simulation data, LPM')
    subplot(224)
        simpleBodephase(sysID.nonPar.OL.frfComp.lpm.raw,'Hz','wrap','b',1);hold on;grid minor
        simpleBodephase(sysID.nonPar.CL.frfComp.lpm.direct.raw ,'Hz','wrap','r',1)
        simpleBodephase(sysID.nonPar.CL.frfComp.lpm.classic.raw,'Hz','wrap','k',1)
        simpleBodephase(sysID.nonPar.CL.frfComp.lpm.coprime.raw,'Hz','wrap','m',1)
        simpleBodephase(G(1,1),'Hz',logspace(-4,0,1000),'wrap','g--')
        xlim([0.0001 0.1])
            legend('OL, LPM','CL, Direct','CL, Indirect classic','CL, Indirect coprime','Data generating model')
    print('./_Figures/nonParam/compensated_FRF_tradLPM_8h','-depsc')


%% Plot the fit and the compensated data for the user to check (seperate plots)
    figure(3103);clf;hold on
    set(gcf,'position',[100,100,500,200])
        plot(timeData,newy);grid minor
        plot(timeData,trainDataFilt.y,'--')
            xlabel('Time [s]')
            ylabel('Temperature [degC]')
            title(['Exponential fit on simulated temperature data'])
            legend('Exponential fit','Simulation data','location','best')
    print('./_Figures/OL/exponentialFit_OL','-depsc')
    figure(3104);clf;hold on
    set(gcf,'position',[100,100,500,200])
        plot(timeData,newy-trainDataFilt.y);grid minor
            xlabel('Time [s]')
            ylabel('Temperature fit error [degC]')
            title(['Exponential fit error on simulated temperature data'])
            legend(['RMSE: ',num2str(rms(newy-trainDataFilt.y)),'degC'])
    print('./_Figures/OL/exponentialFit_diff','-depsc')

%% Monte Carlo of identification methods, both OL and C
clear sysID sysID_MC sys % Make sure a clean workspace is used

compTrans = true; % Transient compensation should be on
newData   = true; % To make a new simulation, should give the same result if the input and noise signal are constant

% Define the input signals for the open loop and closed loop identifications
dist_OL = genMultisine(fs, tMeas*3600*fs, 1, 120,1);
dist_CL = genMultisine(fs, tMeas*3600*fs, 1, 20, 1)+118+25;

% Run the Monte Carlo simulations (20 times)
n_MC = 20;

for i=1:n_MC
    % Define the noise signal, and generate a new one for every simulation
    v = 0.15.*randn(N,1);
    
    % Open loop identification
    dist = dist_OL;
    seperate_OLsim_sysID_V4
    
    % Closed loop identification
    dist = dist_CL;
    seperate_CLsim_sysID_V5_2
    
    % Save the identification for later use
    sysID_MC{i} = sysID;
end

% Posibility to save the data
    % save('./_Data/MonteCarlo/OLCL_nonComp_nMC20_23_Oct_25','sysID_MC','n_MC','G')
    % save('./_Data/MonteCarlo/OLCL_Comp_nMC20_23_Oct_25','sysID_MC','n_MC','G')

%% Define the identified systems differently for easier use
for i=1:n_MC

    sys_OL_fixedOrder5{i} = sysID_MC{i}.par.OL.time.fixedOrder5.prd.filt;
    sys_OL_fixedOrder3{i} = sysID_MC{i}.par.OL.time.fixedOrder.prd.filt;
    sys_OL_fixedOrder1{i} = sysID_MC{i}.par.OL.time.fixedOrder1.prd.filt;
    sys_OL_straight{i}    = sysID_MC{i}.par.OL.time.straight.prd.filt;
    sys_OL_initSys{i}     = sysID_MC{i}.par.OL.time.initSys.prd.filt;
    sys_OL_firstApprox{i} = sysID_MC{i}.par.OL.time.firstApprox.filt;
    
    sys_FRF_OL_trd{i} = sysID_MC{i}.par.OL.freq.trd.prd.filt;
    sys_FRF_OL_lpm{i} = sysID_MC{i}.par.OL.freq.lpm.prd.filt;

    sys_CL_direct{i}   = sysID_MC{i}.par.CL.time.direct.filt;
    sys_CL_coprime{i}  = sysID_MC{i}.par.CL.time.coprime.filt;
    sys_CL_twostage{i} = sysID_MC{i}.par.CL.time.twostage.filt;
end

%% Define the data that can be used for the validation

% Define the data that can be used for the (full data) validation
validData_id_full = sysID_MC{1}.OLdataFilt.full.id;
validData_tVec_full = sysID_MC{1}.OLdata.full.tVec;

% Define the data that can be used for the (normal) validation
validData_id = sysID_MC{1}.OLdataFilt.valid.id;
validData_tVec = sysID_MC{1}.OLdata.valid.tVec;

error_OL_fixedOrder5 = prediction_MC(validData_id,validData_tVec,sys_OL_fixedOrder5);
error_OL_fixedOrder3 = prediction_MC(validData_id,validData_tVec,sys_OL_fixedOrder3);
error_OL_fixedOrder1 = prediction_MC(validData_id,validData_tVec,sys_OL_fixedOrder1);
error_OL_straight    = prediction_MC(validData_id,validData_tVec,sys_OL_straight);
error_OL_initSys     = prediction_MC(validData_id,validData_tVec,sys_OL_initSys);
error_OL_firstApprox = prediction_MC(validData_id,validData_tVec,sys_OL_firstApprox);

error_CL_direct   = prediction_MC(validData_id,validData_tVec,sys_CL_direct  );
error_CL_coprime  = prediction_MC(validData_id,validData_tVec,sys_CL_coprime );
error_CL_twostage = prediction_MC(validData_id,validData_tVec,sys_CL_twostage);

%% Give the RMSE's of the OL predictions
warning off
disp(' ')
disp('Prediction RMSE of Open Loop identification methods, averaged over 20 simulations')
disp(['First Order Approx.: ',num2str(error_OL_firstApprox.avg,4),' [degC]'])
disp(['Initialized System:  ',num2str(error_OL_initSys.avg,4),    ' [degC]'])
disp(['Fixed Order, nx=1:   ',num2str(error_OL_fixedOrder1.avg,4),' [degC]'])
disp(['Fixed Order, nx=3:   ',num2str(error_OL_fixedOrder3.avg,4),' [degC]'])
disp(['Fixed Order, nx=5:   ',num2str(error_OL_fixedOrder5.avg,4),' [degC]'])
disp(['Straight Data:       ',num2str(error_OL_straight.avg,4),   ' [degC]'])

H2_OL_fixedOrder5 = 0;
H2_OL_fixedOrder3 = 0;
H2_OL_fixedOrder1 = 0;
H2_OL_straight    = 0;
H2_OL_initSys     = 0;
H2_OL_firstApprox = 0;
for i=1:n_MC
    H2_OL_fixedOrder5 = H2_OL_fixedOrder5 + 1/n_MC.*norm(sys_OL_firstApprox{i}-G(1,1));
    H2_OL_fixedOrder3 = H2_OL_fixedOrder3 + 1/n_MC.*norm(sys_OL_initSys{i}    -G(1,1));
    H2_OL_fixedOrder1 = H2_OL_fixedOrder1 + 1/n_MC.*norm(sys_OL_fixedOrder1{i}-G(1,1));
    H2_OL_straight    = H2_OL_straight    + 1/n_MC.*norm(sys_OL_fixedOrder3{i}-G(1,1));
    H2_OL_initSys     = H2_OL_initSys     + 1/n_MC.*norm(sys_OL_fixedOrder5{i}-G(1,1));
    H2_OL_firstApprox = H2_OL_firstApprox + 1/n_MC.*norm(sys_OL_straight{i}   -G(1,1));
end

disp(' ')
disp('H2-norm of Open Loop identification methods, averaged over 20 simulations')
disp(['First Order Approx.: ',num2str(H2_OL_firstApprox,4)])
disp(['Initialized System:  ',num2str(H2_OL_initSys    ,4)])
disp(['Fixed Order, nx=1:   ',num2str(H2_OL_fixedOrder1,4)])
disp(['Fixed Order, nx=3:   ',num2str(H2_OL_fixedOrder3,4)])
disp(['Fixed Order, nx=5:   ',num2str(H2_OL_fixedOrder5,4)])
disp(['Straight Data:       ',num2str(H2_OL_straight   ,4)])

gap_OL_fixedOrder5 = 0;
gap_OL_fixedOrder3 = 0;
gap_OL_fixedOrder1 = 0;
gap_OL_straight    = 0;
gap_OL_initSys     = 0;
gap_OL_firstApprox = 0;
for i=1:n_MC
    gap_OL_fixedOrder5 = gap_OL_fixedOrder5 + 1/n_MC.*gapmetric(sys_OL_firstApprox{i},G(1,1));
    gap_OL_fixedOrder3 = gap_OL_fixedOrder3 + 1/n_MC.*gapmetric(sys_OL_initSys{i}    ,G(1,1));
    gap_OL_fixedOrder1 = gap_OL_fixedOrder1 + 1/n_MC.*gapmetric(sys_OL_fixedOrder1{i},G(1,1));
    gap_OL_straight    = gap_OL_straight    + 1/n_MC.*gapmetric(sys_OL_fixedOrder3{i},G(1,1));
    gap_OL_initSys     = gap_OL_initSys     + 1/n_MC.*gapmetric(sys_OL_fixedOrder5{i},G(1,1));
    gap_OL_firstApprox = gap_OL_firstApprox + 1/n_MC.*gapmetric(sys_OL_straight{i}   ,G(1,1));
end

disp(' ')
disp('Gap metric of Open Loop identification methods, averaged over 20 simulations')
disp(['First Order Approx.: ',num2str(gap_OL_firstApprox,4)])
disp(['Initialized System:  ',num2str(gap_OL_initSys    ,4)])
disp(['Fixed Order, nx=1:   ',num2str(gap_OL_fixedOrder1,4)])
disp(['Fixed Order, nx=3:   ',num2str(gap_OL_fixedOrder3,4)])
disp(['Fixed Order, nx=5:   ',num2str(gap_OL_fixedOrder5,4)])
disp(['Straight Data:       ',num2str(gap_OL_straight   ,4)])
            
%% Give the RMSE's of the CL predictions
warning off
disp(' ')
disp('Prediction RMSE of Closed Loop identification methods, averaged over 20 simulations')
disp(['Direct method: ',num2str(error_CL_direct.avg  ,4),' [degC]'])
disp(['Co-prime:      ',num2str(error_CL_coprime.avg ,4),' [degC]'])
disp(['two-stage:     ',num2str(error_CL_twostage.avg,4),' [degC]'])

H2_CL_direct   = 0;
H2_CL_coprime  = 0;
H2_CL_twostage = 0;
for i=1:n_MC
    H2_CL_direct   = H2_CL_direct   + 1/n_MC.*norm(sys_CL_direct{i}  -G(1,1));
    H2_CL_coprime  = H2_CL_coprime  + 1/n_MC.*norm(sys_CL_coprime{i} -G(1,1));
    H2_CL_twostage = H2_CL_twostage + 1/n_MC.*norm(sys_CL_twostage{i}-G(1,1));
end

disp(' ')
disp('H2-norm of Closed Loop identification methods, averaged over 20 simulations')
disp(['First Order Approx.: ',num2str(H2_CL_direct  ,4)])
disp(['Initialized System:  ',num2str(H2_CL_coprime ,4)])
disp(['Fixed Order, nx=1:   ',num2str(H2_CL_twostage,4)])

gap_CL_direct   = 0;
gap_CL_coprime  = 0;
gap_CL_twostage = 0;
for i=1:n_MC
    gap_CL_direct   = gap_CL_direct   + 1/n_MC.*gapmetric(sys_CL_direct{i}  ,G(1,1));
    gap_CL_coprime  = gap_CL_coprime  + 1/n_MC.*gapmetric(sys_CL_coprime{i} ,G(1,1));
    gap_CL_twostage = gap_CL_twostage + 1/n_MC.*gapmetric(sys_CL_twostage{i},G(1,1));
end

disp(' ')
disp('Gap Metric of Closed Loop identification methods, averaged over 20 simulations')
disp(['First Order Approx.: ',num2str(gap_CL_direct  ,4)])
disp(['Initialized System:  ',num2str(gap_CL_coprime ,4)])
disp(['Fixed Order, nx=1:   ',num2str(gap_CL_twostage,4)])

%% Make box plots of prediction errors

% validData_id = sysID_MC{1}.OLdataFilt.valid.id;
% validData_tVec = sysID_MC{1}.OLdata.valid.tVec;
% 
% error_OL_fixedOrder5 = prediction_MC(validData_id,validData_tVec,sys_OL_fixedOrder5);
% error_OL_fixedOrder3 = prediction_MC(validData_id,validData_tVec,sys_OL_fixedOrder3);
% error_OL_fixedOrder1 = prediction_MC(validData_id,validData_tVec,sys_OL_fixedOrder1);
% error_OL_straight    = prediction_MC(validData_id,validData_tVec,sys_OL_straight);
% error_OL_initSys     = prediction_MC(validData_id,validData_tVec,sys_OL_initSys);
% error_OL_firstApprox = prediction_MC(validData_id,validData_tVec,sys_OL_firstApprox);
% 
% error_CL_direct   = prediction_MC(validData_id,validData_tVec,sys_CL_direct  );
% error_CL_coprime  = prediction_MC(validData_id,validData_tVec,sys_CL_coprime );
% error_CL_twostage = prediction_MC(validData_id,validData_tVec,sys_CL_twostage);

error_all_OL = [error_OL_fixedOrder5.all;...
                error_OL_fixedOrder3.all;...
                error_OL_fixedOrder1.all;...
                error_OL_straight.all;...
                error_OL_initSys.all;...
                error_OL_firstApprox.all];
errorNames_OL = {'Fixed Order, nx=5';...
                 'Fixed Order, nx=3';...
                 'Fixed Order, nx=1';...
                 'Straight Data';...
                 'Initialized System';...
                 'First Order Approx.'};
error_all_OL_zoomed1 = [error_OL_fixedOrder5.all;...
                       error_OL_fixedOrder3.all];
errorNames_OL_zoomed1 = {'Fixed Order, nx=5';...
                        'Fixed Order, nx=3'};
error_all_OL_zoomed2 = [error_OL_initSys.all];
errorNames_OL_zoomed2 = {'Initialized System'};
% error_all_OL_zoomed = [error_OL_fixedOrder5.all;...
%                        error_OL_fixedOrder3.all];
% errorNames_OL_zoomed = {'Fixed Order, nx=5';...
%                         'Fixed Order, nx=3'};

error_all_CL = [error_CL_direct.all  ;...
                error_CL_coprime.all ;...
                error_CL_twostage.all];
errorNames_CL = {'Closed Loop, Direct'   ;...
                 'Closed Loop, Coprime'  ;...
                 'Closed Loop, Two-stage'};
error_all_CL_zoomed = [error_CL_direct.all];
errorNames_CL_zoomed = {'Closed Loop, Direct'};

figure(3201);clf
set(gcf,'position',[100 100 600 250])
    boxplot(error_all_OL',errorNames_OL);grid minor
        ylabel('Prediction RMSE [K]')
        title('Prediction RMSE distribution, all OL methods')
    print('.\_Figures\OL\parTime/boxplots_MC_predRMSE_simOLmethods_8h','-depsc')

figure(3202);clf
set(gcf,'position',[100 100 600 250])
    boxplot(error_all_CL',errorNames_CL);grid minor
        ylabel('Prediction RMSE [K]')
        title('Prediction RMSE distribution, all CL methods')
    print('.\_Figures\CL\parTime/boxplots_MC_predRMSE_simCLmethods_8h','-depsc')


% figure();clf
%     subplot(121)
%         boxplot(error_all_OL_zoomed1',errorNames_OL_zoomed1);grid minor
%             title('Prediction RMSE distribution, selected OL methods')
%     subplot(122)
%         boxplot(error_all_OL_zoomed2',errorNames_OL_zoomed2);grid minor
%             title('Prediction RMSE distribution, selected OL methods')
% 
% figure()
%     boxplot(error_all_CL_zoomed',errorNames_CL_zoomed);grid minor
%         title('Prediction RMSE distribution, selected CL methods')

%% Visualize the Monte Carlo sim, for parametric models

figure(3212);clf
sgtitle('Monte Carlo on OL identification, 20 Simulations')
    subplot(431)
        simpleBodemag_MC(sys_OL_fixedOrder5,'Hz',bodeRange,'b')
        simpleBodemag(G(1,1),'Hz',bodeRange,'k--',lw)
            title('Fixed Order, nx=5')
    subplot(434)
        simpleBodephase_MC(sys_OL_fixedOrder5,'Hz','wrap',bodeRange,'b')
        simpleBodephase(G(1,1),'Hz','wrap',bodeRange,'k--',lw)
    subplot(432)
        simpleBodemag_MC(sys_OL_fixedOrder3,'Hz',bodeRange,'b')
        simpleBodemag(G(1,1),'Hz',bodeRange,'k--',lw)
            title('Fixed Order, nx=3')
    subplot(435)
        simpleBodephase_MC(sys_OL_fixedOrder3,'Hz','wrap',bodeRange,'b')
        simpleBodephase(G(1,1),'Hz','wrap',bodeRange,'k--',lw)
    subplot(433)
        simpleBodemag_MC(sys_OL_fixedOrder1,'Hz',bodeRange,'b')
        simpleBodemag(G(1,1),'Hz',bodeRange,'k--',lw)
            title('Fixed Order, nx=1')
    subplot(436)
        simpleBodephase_MC(sys_OL_fixedOrder1,'Hz','wrap',bodeRange,'b')
        simpleBodephase(G(1,1),'Hz','wrap',bodeRange,'k--',lw)
    subplot(437)
        simpleBodemag_MC(sys_OL_straight,'Hz',bodeRange,'b')
        simpleBodemag(G(1,1),'Hz',bodeRange,'k--',lw)
            title('Straight')
    subplot(4,3,10)
        simpleBodephase_MC(sys_OL_straight,'Hz','wrap',bodeRange,'b')
        simpleBodephase(G(1,1),'Hz','wrap',bodeRange,'k--',lw)
    subplot(438)
        simpleBodemag_MC(sys_OL_initSys,'Hz',bodeRange,'b')
        simpleBodemag(G(1,1),'Hz',bodeRange,'k--',lw)
            title('Initialized system')
    subplot(4,3,11)
        simpleBodephase_MC(sys_OL_initSys,'Hz','wrap',bodeRange,'b')
        simpleBodephase(G(1,1),'Hz','wrap',bodeRange,'k--',lw)
    subplot(439)
        simpleBodemag_MC(sys_OL_firstApprox,'Hz',bodeRange,'b')
        simpleBodemag(G(1,1),'Hz',bodeRange,'k--',lw)
            title('First Approximation')
    subplot(4,3,12)
        simpleBodephase_MC(sys_OL_firstApprox,'Hz','wrap',bodeRange,'b')
        simpleBodephase(G(1,1),'Hz','wrap',bodeRange,'k--',lw)
    print('./_Figures/CL/parTime/variance_MC_bodes_simCLmethods_8h','-depsc')

figure(3211);clf
sgtitle('Monte Carlo on CL identification, 20 Simulations')
    subplot(231)
        simpleBodemag_MC(sys_CL_direct,'Hz',bodeRange,'b')
        simpleBodemag(G(1,1),'Hz',bodeRange,'k--',lw)
            title('Direct')
    subplot(234)
        simpleBodephase_MC(sys_CL_direct,'Hz','wrap',bodeRange,'b')
        simpleBodephase(G(1,1),'Hz','wrap',bodeRange,'k--',lw)
    subplot(232)
        simpleBodemag_MC(sys_CL_coprime,'Hz',bodeRange,'b')
        simpleBodemag(G(1,1),'Hz',bodeRange,'k--',lw)
            title('Co-prime')
    subplot(235)
        simpleBodephase_MC(sys_CL_coprime,'Hz','wrap',bodeRange,'b')
        simpleBodephase(G(1,1),'Hz','wrap',bodeRange,'k--',lw)
    subplot(233)
        simpleBodemag_MC(sys_CL_twostage,'Hz',bodeRange,'b')
        simpleBodemag(G(1,1),'Hz',bodeRange,'k--',lw)
            title('Two-stage')
    subplot(236)
        simpleBodephase_MC(sys_CL_twostage,'Hz','wrap',bodeRange,'b')
        simpleBodephase(G(1,1),'Hz','wrap',bodeRange,'k--',lw)
    print('./_Figures/CL/parTime/variance_MC_bodes_simCLmethods_8h','-depsc')

%% Visualize the Monte Carlo sim, for non-parametric models
for i=1:n_MC
    sys1{i} = sysID_MC{i}.nonPar.OL.frfComp.trd.filt;
    sys2{i} = sysID_MC{i}.nonPar.OL.frfComp.lpm.filt;
    sys3{i} = sysID_MC{i}.nonPar.CL.frfComp.trd.direct.filt;
    sys4{i} = sysID_MC{i}.nonPar.CL.frfComp.lpm.direct.filt;
end

figure(3213);clf
    subplot(211)
        simpleBodemag_MC(sys1,'Hz','b');grid minor
        simpleBodemag_MC(sys2,'Hz','r')
        simpleBodemag(G(1,1),'Hz','k--',lw)
        xlim([1e-4 1e-1])
        ylim([-130 0])
    subplot(212)
        simpleBodephase_MC(sys1,'Hz','wrap','b');grid minor
        simpleBodephase_MC(sys2,'Hz','wrap','r')
        simpleBodephase(G(1,1),'Hz','wrap','k--',lw)
        xlim([1e-4 1e-1])
    print('./_Figures/OL/nonPar/variance_MC_FRF_simOLmethods_8h','-depsc')

figure(3214);clf
    subplot(211)
        simpleBodemag_MC(sys3,'Hz','b');grid minor;hold on
        simpleBodemag_MC(sys4,'Hz','r')
        simpleBodemag(G(1,1),'Hz','k--',lw)
        xlim([1e-4 1e-1])
        ylim([-130 0])
    subplot(212)
        simpleBodephase_MC(sys3,'Hz','wrap','b');grid minor;hold on
        simpleBodephase_MC(sys4,'Hz','wrap','r')
        simpleBodephase(G(1,1),'Hz','wrap','k--',lw)
        xlim([1e-4 1e-1])
    print('./_Figures/CL/nonPar/variance_MC_FRF_simCLmethods_8h','-depsc')

%% Monte Carlo of open loop identification, different data lengths
clear sysID sysID_MC sys
compTrans = true;
newData   = true;

times = 8:2:24;

dist_OL = genMultisine(fs, times(end)*3600*fs, 1, 120,1);
dist_CL = genMultisine(fs, times(end)*3600*fs, 1, 20, 1)+118+25;

n_MC = 10;
for i=1:n_MC
    v_full = 0.15.*randn((times(end)*3600*fs),1);
    for j=1:length(times)
        v = v_full(1:times(j)*3600*fs);
        tMeas = times(j);

        dist = dist_OL(1:(times(j)*3600*fs));
        seperate_OLsim_sysID_V4
        
        dist = dist_CL(1:(times(j)*3600*fs));
        seperate_CLsim_sysID_V5_2
        
        sysID_MC{i,j} = sysID;
    end
end

% Posibility to save data again
    % save('./_Data/MonteCarlo/increaseTime_OLCL_Comp_nMC3_24_Oct_25','sysID_MC','n_MC','G')

%% Define the used systems
for i=1:n_MC
    for j=1:length(times)
    sys_time_OL_fixedOrder5{i,j} = sysID_MC{i,j}.par.OL.time.fixedOrder5.prd.filt;
    sys_time_OL_fixedOrder3{i,j} = sysID_MC{i,j}.par.OL.time.fixedOrder.prd.filt;
    sys_time_OL_fixedOrder1{i,j} = sysID_MC{i,j}.par.OL.time.fixedOrder1.prd.filt;
    sys_time_OL_straight{i,j}    = sysID_MC{i,j}.par.OL.time.straight.prd.filt;
    sys_time_OL_initSys{i,j}     = sysID_MC{i,j}.par.OL.time.initSys.prd.filt;
    sys_time_OL_firstApprox{i,j} = sysID_MC{i,j}.par.OL.time.firstApprox.filt;
    
    sys_time_FRF_OL_trd{i} = sysID_MC{i}.par.OL.freq.trd.prd.filt;
    sys_time_FRF_OL_lpm{i} = sysID_MC{i}.par.OL.freq.lpm.prd.filt;

    sys_time_CL_direct{i,j}   = sysID_MC{i,j}.par.CL.time.direct.filt;
    sys_time_CL_coprime{i,j}  = sysID_MC{i,j}.par.CL.time.coprime.filt;
    sys_time_CL_twostage{i,j} = sysID_MC{i,j}.par.CL.time.twostage.filt;
    end
end

for j=1:length(times)
    error_time_OL_fixedOrder5{j} = prediction_MC(validData_id,validData_tVec,sys_time_OL_fixedOrder5(:,j)');
    error_time_OL_fixedOrder3{j} = prediction_MC(validData_id,validData_tVec,sys_time_OL_fixedOrder3(:,j)');
    error_time_OL_fixedOrder1{j} = prediction_MC(validData_id,validData_tVec,sys_time_OL_fixedOrder1(:,j)');
    error_time_OL_straight{j}    = prediction_MC(validData_id,validData_tVec,sys_time_OL_straight(:,j)');
    error_time_OL_initSys{j}     = prediction_MC(validData_id,validData_tVec,sys_time_OL_initSys(:,j)');
    error_time_OL_firstApprox{j} = prediction_MC(validData_id,validData_tVec,sys_time_OL_firstApprox(:,j)');
    
    error_time_CL_direct{j}   = prediction_MC(validData_id,validData_tVec,sys_time_CL_direct(:,j)'  );
    error_time_CL_coprime{j}  = prediction_MC(validData_id,validData_tVec,sys_time_CL_coprime(:,j)' );
    error_time_CL_twostage{j} = prediction_MC(validData_id,validData_tVec,sys_time_CL_twostage(:,j)');
end

%% Plot average RMSE over measurement time

for j=1:length(times)
    RMSE_time_OL_fixedOrder5(j) = error_time_OL_fixedOrder5{j}.avg.*1000;
    RMSE_time_OL_fixedOrder3(j) = error_time_OL_fixedOrder3{j}.avg.*1000;
    RMSE_time_OL_fixedOrder1(j) = error_time_OL_fixedOrder1{j}.avg.*1000;
    RMSE_time_OL_straight(j)    = error_time_OL_straight{j}.avg.*1000;
    RMSE_time_OL_initSys(j)     = error_time_OL_initSys{j}.avg.*1000;
    RMSE_time_OL_firstApprox(j) = error_time_OL_firstApprox{j}.avg.*1000;

    RMSEall_time_OL_fixedOrder5(j,:) = error_time_OL_fixedOrder5{j}.all;
    RMSEall_time_OL_fixedOrder3(j,:) = error_time_OL_fixedOrder3{j}.all;
    RMSEall_time_OL_fixedOrder1(j,:) = error_time_OL_fixedOrder1{j}.all;
    RMSEall_time_OL_straight(j,:)    = error_time_OL_straight{j}.all;
    RMSEall_time_OL_initSys(j,:)     = error_time_OL_initSys{j}.all;
    RMSEall_time_OL_firstApprox(j,:) = error_time_OL_firstApprox{j}.all;

    RMSE_time_CL_direct(j)   = error_time_CL_direct{j}.avg.*1000;
    RMSE_time_CL_coprime(j)  = error_time_CL_coprime{j}.avg.*1000;
    RMSE_time_CL_twostage(j) = error_time_CL_twostage{j}.avg.*1000;

    RMSEall_time_CL_direct(j,:)   = error_time_CL_direct{j}.all;
    RMSEall_time_CL_coprime(j,:)  = error_time_CL_coprime{j}.all;
    RMSEall_time_CL_twostage(j,:) = error_time_CL_twostage{j}.all;
end

max_RMSEall_time_OL_fixedOrder5 = max(RMSEall_time_OL_fixedOrder5')'.*1000;
max_RMSEall_time_OL_fixedOrder3 = max(RMSEall_time_OL_fixedOrder3')'.*1000;
max_RMSEall_time_OL_fixedOrder1 = max(RMSEall_time_OL_fixedOrder1')'.*1000;
max_RMSEall_time_OL_straight    = max(RMSEall_time_OL_straight'   )'.*1000;
max_RMSEall_time_OL_initSys     = max(RMSEall_time_OL_initSys'    )'.*1000;
max_RMSEall_time_OL_firstApprox = max(RMSEall_time_OL_firstApprox')'.*1000;

min_RMSEall_time_OL_fixedOrder5 = min(RMSEall_time_OL_fixedOrder5')'.*1000;
min_RMSEall_time_OL_fixedOrder3 = min(RMSEall_time_OL_fixedOrder3')'.*1000;
min_RMSEall_time_OL_fixedOrder1 = min(RMSEall_time_OL_fixedOrder1')'.*1000;
min_RMSEall_time_OL_straight    = min(RMSEall_time_OL_straight'   )'.*1000;
min_RMSEall_time_OL_initSys     = min(RMSEall_time_OL_initSys'    )'.*1000;
min_RMSEall_time_OL_firstApprox = min(RMSEall_time_OL_firstApprox')'.*1000;


max_RMSEall_time_CL_direct   = max(RMSEall_time_CL_direct'  )'.*1000;
max_RMSEall_time_CL_coprime  = max(RMSEall_time_CL_coprime' )'.*1000;
max_RMSEall_time_CL_twostage = max(RMSEall_time_CL_twostage')'.*1000;

min_RMSEall_time_CL_direct   = min(RMSEall_time_CL_direct'  )'.*1000;
min_RMSEall_time_CL_coprime  = min(RMSEall_time_CL_coprime' )'.*1000;
min_RMSEall_time_CL_twostage = min(RMSEall_time_CL_twostage')'.*1000;


figure(3221);clf
set(gcf,'position',[100,100,1300 450])
    subplot(235);hold on;grid minor
        fill([times flip(times)],[min_RMSEall_time_OL_fixedOrder5;flip(max_RMSEall_time_OL_fixedOrder5)],[0.8 0.8 1], 'EdgeColor', 'none');
        plot(times,RMSE_time_OL_fixedOrder5);
            xlabel('Available time [h]')
            ylabel('Prediction RMSE [mK]')
            title('Average and spread of RMSE, Fixed Order, nx=5')
    subplot(234);hold on;grid minor
        fill([times flip(times)],[min_RMSEall_time_OL_fixedOrder3;flip(max_RMSEall_time_OL_fixedOrder3)],[0.8 0.8 1], 'EdgeColor', 'none');
        plot(times,RMSE_time_OL_fixedOrder3);
            xlabel('Available time [h]')
            ylabel('Prediction RMSE [mK]')
            title('Average and spread of RMSE, Fixed Order, nx=3')
    subplot(233);hold on;grid minor
        fill([times flip(times)],[min_RMSEall_time_OL_fixedOrder1;flip(max_RMSEall_time_OL_fixedOrder1)],[0.8 0.8 1], 'EdgeColor', 'none');
        plot(times,RMSE_time_OL_fixedOrder1);
            xlabel('Available time [h]')
            ylabel('Prediction RMSE [mK]')
            title('Average and spread of RMSE, Fixed Order, nx=1')
    subplot(236);hold on;grid minor
        fill([times flip(times)],[min_RMSEall_time_OL_straight;flip(max_RMSEall_time_OL_straight)],[0.8 0.8 1], 'EdgeColor', 'none');
        plot(times,RMSE_time_OL_straight);
            xlabel('Available time [h]')
            ylabel('Prediction RMSE [mK]')
            title('Average and spread of RMSE, Straight data')
    subplot(232);hold on;grid minor
        fill([times flip(times)],[min_RMSEall_time_OL_initSys;flip(max_RMSEall_time_OL_initSys)],[0.8 0.8 1], 'EdgeColor', 'none');
        plot(times,RMSE_time_OL_initSys);
            xlabel('Available time [h]')
            ylabel('Prediction RMSE [mK]')
            title('Average and spread of RMSE, Initialized system')
    subplot(231);hold on;grid minor
        fill([times flip(times)],[min_RMSEall_time_OL_firstApprox;flip(max_RMSEall_time_OL_firstApprox)],[0.8 0.8 1], 'EdgeColor', 'none');
        plot(times,RMSE_time_OL_firstApprox);
            xlabel('Available time [h]')
            ylabel('Prediction RMSE [mK]')
            title('Average and spread of RMSE, First approx.')
    print('./_Figures/OL/parTime/timeVariance_MC_bodes_simOLmethods_upTo24h','-depsc')

figure(3222);clf
set(gcf,'position',[100,100,1300 300])
    subplot(131);hold on;grid minor
        fill([times flip(times)],[min_RMSEall_time_CL_direct;flip(max_RMSEall_time_CL_direct)],[0.8 0.8 1], 'EdgeColor', 'none');
        plot(times,RMSE_time_CL_direct);
            xlabel('Available time [h]')
            ylabel('Prediction RMSE [mK]')
            title('Average and spread of RMSE, Direct method')
    subplot(132);hold on;grid minor
        fill([times flip(times)],[min_RMSEall_time_CL_coprime;flip(max_RMSEall_time_CL_coprime)],[0.8 0.8 1], 'EdgeColor', 'none');
        plot(times,RMSE_time_CL_coprime);
            xlabel('Available time [h]')
            ylabel('Prediction RMSE [mK]')
            title('Average and spread of RMSE, Co-prime method')
    subplot(133);hold on;grid minor
        fill([times flip(times)],[min_RMSEall_time_CL_twostage;flip(max_RMSEall_time_CL_twostage)],[0.8 0.8 1], 'EdgeColor', 'none');
        plot(times,RMSE_time_CL_twostage);
            xlabel('Available time [h]')
            ylabel('Prediction RMSE [mK]')
            title('Average and spread of RMSE, Two-stage method')
    print('./_Figures/CL/parTime/timeVariance_MC_bodes_simCLmethods_upTo24h','-depsc')





























%%
%% Might still be usefull sometime, but are not in the report, so not needed right now
%%
%% Compare the parametric identification of OL data

figure(301);clf
set(gcf,'position',[300 100 700 450])
    subplot(211)
        simpleBodemag(sysID.par.OL.time.firstApprox.filt    ,'Hz',1,bodeRange,'b');hold on;grid minor
        simpleBodemag(sysID.par.OL.time.initSys.prd.filt    ,'Hz',1,bodeRange,'r')
        simpleBodemag(sysID.par.OL.time.fixedOrder1.prd.filt,'Hz',1,bodeRange,'k')
        simpleBodemag(sysID.par.OL.time.fixedOrder5.prd.filt,'Hz',1,bodeRange,'m') 
        simpleBodemag(sysID.par.OL.time.straight.prd.filt   ,'Hz',1,bodeRange,'c')
        simpleBodemag(G(1,1)-sysID.par.OL.time.firstApprox.filt    ,'Hz',1,bodeRange,'b-.')
        simpleBodemag(G(1,1)-sysID.par.OL.time.initSys.prd.filt    ,'Hz',1,bodeRange,'r-.')
        simpleBodemag(G(1,1)-sysID.par.OL.time.fixedOrder1.prd.filt,'Hz',1,bodeRange,'k-.')
        simpleBodemag(G(1,1)-sysID.par.OL.time.fixedOrder5.prd.filt,'Hz',1,bodeRange,'m-.') 
        simpleBodemag(G(1,1)-sysID.par.OL.time.straight.prd.filt   ,'Hz',1,bodeRange,'c-.')
        simpleBodemag(G(1,1)                                ,'Hz',1,bodeRange,'g--')
            xlim([bodeRange(1) bodeRange(end)])
            title('Bode plots of a parametric identification methods using filtered open loop data')
    subplot(212)
        simpleBodephase(sysID.par.OL.time.firstApprox.filt    ,'Hz',1,'wrap',bodeRange,'b');hold on;grid minor
        simpleBodephase(sysID.par.OL.time.initSys.prd.filt    ,'Hz',1,'wrap',bodeRange,'r')
        simpleBodephase(sysID.par.OL.time.fixedOrder1.prd.filt,'Hz',1,'wrap',bodeRange,'k')
        simpleBodephase(sysID.par.OL.time.fixedOrder5.prd.filt,'Hz',1,'wrap',bodeRange,'m')
        simpleBodephase(sysID.par.OL.time.straight.prd.filt   ,'Hz',1,'wrap',bodeRange,'c')
        simpleBodephase(G(1,1)                                ,'Hz',1,'wrap',bodeRange,'g--')
            xlim([bodeRange(1) bodeRange(end)])
            legend('First order approx.','Initialized system','Fixed order, nx=1','Fixed order, nx=5','Straight data','Data generating model')
print('./_Figures/OL/parTime/comparison_simOLmethods_8h','-depsc')

%% Calculate the gap metric of the identified systems wrt the data generating model
[sysID.gap.OL.time.firstApprox.filt    ,~] = gapmetric(G(1,1),sysID.par.OL.time.firstApprox.filt    );
[sysID.gap.OL.time.initSys.prd.filt    ,~] = gapmetric(G(1,1),sysID.par.OL.time.initSys.prd.filt    );
[sysID.gap.OL.time.fixedOrder1.prd.filt,~] = gapmetric(G(1,1),sysID.par.OL.time.fixedOrder1.prd.filt);
[sysID.gap.OL.time.fixedOrder5.prd.filt,~] = gapmetric(G(1,1),sysID.par.OL.time.fixedOrder5.prd.filt);
[sysID.gap.OL.time.straight.prd.filt   ,~] = gapmetric(G(1,1),sysID.par.OL.time.straight.prd.filt   );

%%
disp(' ')
disp('Gap metric of open loop identification compared to data generating model: ')
disp(['First order approximation: ',num2str(sysID.gap.OL.time.firstApprox.filt    )])
disp(['Initialized system:        ',num2str(sysID.gap.OL.time.initSys.prd.filt    )])
disp(['Fixed order, nx=1:         ',num2str(sysID.gap.OL.time.fixedOrder1.prd.filt)])
disp(['Fixed order, nx=5:         ',num2str(sysID.gap.OL.time.fixedOrder5.prd.filt)])
disp(['Straight data:             ',num2str(sysID.gap.OL.time.straight.prd.filt   )])
disp(' ')
disp('H2 norm of open loop identification compared to data generating model: ')
disp(['First order approximation: ',num2str(norm(G(1,1)-sysID.par.OL.time.firstApprox.filt    ))])
disp(['Initialized system:        ',num2str(norm(G(1,1)-sysID.par.OL.time.initSys.prd.filt    ))])
disp(['Fixed order, nx=1:         ',num2str(norm(G(1,1)-sysID.par.OL.time.fixedOrder1.prd.filt))])
disp(['Fixed order, nx=5:         ',num2str(norm(G(1,1)-sysID.par.OL.time.fixedOrder5.prd.filt))])
disp(['Straight data:             ',num2str(norm(G(1,1)-sysID.par.OL.time.straight.prd.filt   ))])

%% Calculate the RMSE of 1-step ahead prediction of the identified systems 
[sysID.pred.OL.time.firstApprox.filt    ] = compare(sysID.OLdataFilt.valid.id,sysID.par.OL.time.firstApprox.filt    ,1);
[sysID.pred.OL.time.initSys.prd.filt    ] = compare(sysID.OLdataFilt.valid.id,sysID.par.OL.time.initSys.prd.filt    ,1);
[sysID.pred.OL.time.fixedOrder1.prd.filt] = compare(sysID.OLdataFilt.valid.id,sysID.par.OL.time.fixedOrder1.prd.filt,1);
[sysID.pred.OL.time.fixedOrder5.prd.filt] = compare(sysID.OLdataFilt.valid.id,sysID.par.OL.time.fixedOrder5.prd.filt,1);
[sysID.pred.OL.time.straight.prd.filt   ] = compare(sysID.OLdataFilt.valid.id,sysID.par.OL.time.straight.prd.filt   ,1);

[sysID.predErr.OL.time.firstApprox.filt    ] = sysID.OLdataFilt.valid.y-sysID.pred.OL.time.firstApprox.filt.y    ;
[sysID.predErr.OL.time.initSys.prd.filt    ] = sysID.OLdataFilt.valid.y-sysID.pred.OL.time.initSys.prd.filt.y    ;
[sysID.predErr.OL.time.fixedOrder1.prd.filt] = sysID.OLdataFilt.valid.y-sysID.pred.OL.time.fixedOrder1.prd.filt.y;
[sysID.predErr.OL.time.fixedOrder5.prd.filt] = sysID.OLdataFilt.valid.y-sysID.pred.OL.time.fixedOrder5.prd.filt.y;
[sysID.predErr.OL.time.straight.prd.filt   ] = sysID.OLdataFilt.valid.y-sysID.pred.OL.time.straight.prd.filt.y   ;

figure(302);clf
    ax1 = subplot(211);
        plot(sysID.OLdata.valid.tVec,sysID.pred.OL.time.firstApprox.filt.y    ,'b',LineWidth=lw);hold on;grid minor
        plot(sysID.OLdata.valid.tVec,sysID.pred.OL.time.initSys.prd.filt.y    ,'r',LineWidth=lw)
        plot(sysID.OLdata.valid.tVec,sysID.pred.OL.time.fixedOrder1.prd.filt.y,'k',LineWidth=lw)
        plot(sysID.OLdata.valid.tVec,sysID.pred.OL.time.fixedOrder5.prd.filt.y,'m',LineWidth=lw)
        plot(sysID.OLdata.valid.tVec,sysID.pred.OL.time.straight.prd.filt.y   ,'c',LineWidth=lw)
        plot(sysID.OLdata.valid.tVec,sysID.OLdataFilt.valid.y               ,'g--',LineWidth=lw)
            xlim([sysID.OLdata.valid.tVec(1) sysID.OLdata.valid.tVec(end)])
            xlabel('Time [s]')
            ylabel('Temperature [degC]')
            title('1-Step ahead prediction of open loop identified models')
            legend('First order approx.','Initialized system','Fixed order, nx=1','Fixed order, nx=5','Straight data','Validation data')
    ax2 = subplot(212);
        plot(sysID.OLdata.valid.tVec,sysID.predErr.OL.time.firstApprox.filt    ,'b',LineWidth=lw);hold on;grid minor
        plot(sysID.OLdata.valid.tVec,sysID.predErr.OL.time.initSys.prd.filt    ,'r',LineWidth=lw)
        plot(sysID.OLdata.valid.tVec,sysID.predErr.OL.time.fixedOrder1.prd.filt,'k',LineWidth=lw)
        plot(sysID.OLdata.valid.tVec,sysID.predErr.OL.time.fixedOrder5.prd.filt,'m',LineWidth=lw)
        plot(sysID.OLdata.valid.tVec,sysID.predErr.OL.time.straight.prd.filt   ,'c',LineWidth=lw)
            xlim([sysID.OLdata.valid.tVec(1) sysID.OLdata.valid.tVec(end)])
            title('1-Step ahead prediction errors of open loop identified models')
            xlabel('Time [s]')
            ylabel('Temperature error [degC]')
    linkaxes([ax1 ax2],'x')
    print('./_Figures/OL/parTime/comparison_predictions_simOLmethods_8h','-depsc')

disp(' ')
disp('RMS prediction error of open loop identification: ')
disp(['First order approximation: ',num2str(rms(sysID.predErr.OL.time.firstApprox.filt    )),'[degC]'])
disp(['Initialized system:        ',num2str(rms(sysID.predErr.OL.time.initSys.prd.filt    )),'[degC]'])
disp(['Fixed order, nx=1:         ',num2str(rms(sysID.predErr.OL.time.fixedOrder1.prd.filt)),'[degC]'])
disp(['Fixed order, nx=5:         ',num2str(rms(sysID.predErr.OL.time.fixedOrder5.prd.filt)),'[degC]'])
disp(['Straight data:             ',num2str(rms(sysID.predErr.OL.time.straight.prd.filt   )),'[degC]'])



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Closed Loop %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compare the parametric identification of CL data

figure(401);clf
set(gcf,'position',[300 100 700 450])
    subplot(211)
        simpleBodemag(sysID.par.CL.time.direct.filt  ,'Hz',1,bodeRange,'b');hold on;grid minor
        simpleBodemag(sysID.par.CL.time.classic.filt ,'Hz',1,bodeRange,'r')
        simpleBodemag(sysID.par.CL.time.coprime.filt ,'Hz',1,bodeRange,'k')
        simpleBodemag(sysID.par.CL.time.twostage.filt,'Hz',1,bodeRange,'m') 
        simpleBodemag(G(1,1)-sysID.par.CL.time.direct.filt  ,'Hz',1,bodeRange,'b-.')
        simpleBodemag(G(1,1)-sysID.par.CL.time.classic.filt ,'Hz',1,bodeRange,'r-.')
        simpleBodemag(G(1,1)-sysID.par.CL.time.coprime.filt ,'Hz',1,bodeRange,'k-.')
        simpleBodemag(G(1,1)-sysID.par.CL.time.twostage.filt,'Hz',1,bodeRange,'m-.') 
        simpleBodemag(G(1,1)                                ,'Hz',1,bodeRange,'g--')
            xlim([bodeRange(1) bodeRange(end)])
            title('Bode plots of a parametric identification methods using filtered open loop data')
    subplot(212)
        simpleBodephase(sysID.par.CL.time.direct.filt  ,'Hz',1,'wrap',bodeRange,'b');hold on;grid minor
        simpleBodephase(sysID.par.CL.time.classic.filt ,'Hz',1,'wrap',bodeRange,'r')
        simpleBodephase(sysID.par.CL.time.coprime.filt ,'Hz',1,'wrap',bodeRange,'k')
        simpleBodephase(sysID.par.CL.time.twostage.filt,'Hz',1,'wrap',bodeRange,'m') 
        simpleBodephase(G(1,1)                                ,'Hz',1,'wrap',bodeRange,'g--')
            xlim([bodeRange(1) bodeRange(end)])
            legend('Direct method','Indirect, classic','Indirect, co-prime','Indirect, two-stage','Data generating model')
print('./_Figures/OL/parTime/comparison_simCLmethods_8h','-depsc')

%% Calculate the gap metric of the identified systems wrt the data generating model
[sysID.gap.CL.time.direct.filt  ,~] = gapmetric(G(1,1),sysID.par.CL.time.direct.filt  );
[sysID.gap.CL.time.coprime.filt ,~] = gapmetric(G(1,1),sysID.par.CL.time.coprime.filt );
[sysID.gap.CL.time.twostage.filt,~] = gapmetric(G(1,1),sysID.par.CL.time.twostage.filt);

%%
disp(' ')
disp('Gap metric of closed loop identification compared to data generating model: ')
disp(['Direct method:       ',num2str(sysID.gap.CL.time.direct.filt  )])
disp(['Indirect, co-prime:  ',num2str(sysID.gap.CL.time.coprime.filt )])
disp(['Indirect, two-stage: ',num2str(sysID.gap.CL.time.twostage.filt)])
disp(' ')
disp('H2 norm of closed loop identification compared to data generating model: ')
disp(['Direct method:       ',num2str(norm(G(1,1)-sysID.gap.CL.time.direct.filt  ))])
disp(['Indirect, co-prime:  ',num2str(norm(G(1,1)-sysID.gap.CL.time.coprime.filt ))])
disp(['Indirect, two-stage: ',num2str(norm(G(1,1)-sysID.gap.CL.time.twostage.filt))])

%% Calculate the RMSE of 1-step ahead prediction of the identified systems
[sysID.pred.CL.time.direct.filt  ] = compare(sysID.CLdataFilt.valid.direct,sysID.par.CL.time.direct.filt  ,1);
[sysID.pred.CL.time.coprime.filt ] = compare(sysID.CLdataFilt.valid.direct,sysID.par.CL.time.coprime.filt ,1);
[sysID.pred.CL.time.twostage.filt] = compare(sysID.CLdataFilt.valid.direct,sysID.par.CL.time.twostage.filt,1);

[sysID.predErr.CL.time.direct.filt  ] = sysID.CLdataFilt.valid.y-sysID.pred.CL.time.direct.filt.y  ;
[sysID.predErr.CL.time.coprime.filt ] = sysID.CLdataFilt.valid.y-sysID.pred.CL.time.coprime.filt.y ;
[sysID.predErr.CL.time.twostage.filt] = sysID.CLdataFilt.valid.y-sysID.pred.CL.time.twostage.filt.y;

figure(402);clf
    ax1 = subplot(211);
        plot(sysID.CLdata.valid.tVec,sysID.pred.CL.time.direct.filt.y  ,'b',LineWidth=lw);hold on;grid minor
        plot(sysID.CLdata.valid.tVec,sysID.pred.CL.time.coprime.filt.y ,'k',LineWidth=lw)
        plot(sysID.CLdata.valid.tVec,sysID.pred.CL.time.twostage.filt.y,'m',LineWidth=lw)
        plot(sysID.CLdata.valid.tVec,sysID.CLdataFilt.valid.y               ,'g--',LineWidth=lw)
            xlim([sysID.CLdata.valid.tVec(1) sysID.CLdata.valid.tVec(end)])
            xlabel('Time [s]')
            ylabel('Temperature [degC]')
            title('1-Step ahead prediction of closed loop identified models')
            legend('Direct method','Indirect, co-prime','Indirect, two-stage','Validation data')
    ax2 = subplot(212);
        plot(sysID.CLdata.valid.tVec,sysID.predErr.CL.time.direct.filt  ,'b',LineWidth=lw);hold on;grid minor
        plot(sysID.CLdata.valid.tVec,sysID.predErr.CL.time.coprime.filt ,'k',LineWidth=lw)
        plot(sysID.CLdata.valid.tVec,sysID.predErr.CL.time.twostage.filt,'m',LineWidth=lw)
            xlim([sysID.CLdata.valid.tVec(1) sysID.CLdata.valid.tVec(end)])
            title('1-Step ahead prediction errors of closed loop identified models')
            xlabel('Time [s]')
            ylabel('Temperature error [degC]')
    linkaxes([ax1 ax2],'x')
    print('./_Figures/OL/parTime/comparison_predictions_simCLmethods_8h','-depsc')

disp(' ')
disp('RMS prediction error of closed loop identification: ')
disp(['Direct method:       ',num2str(rms(sysID.predErr.CL.time.direct.filt  )),'[degC]'])
disp(['Indirect, co-prime:  ',num2str(rms(sysID.predErr.CL.time.coprime.filt )),'[degC]'])
disp(['Indirect, two-stage: ',num2str(rms(sysID.predErr.CL.time.twostage.filt)),'[degC]'])

%% Neural network training on residuals using open loop models

% Model = sysID.par.OL.time.firstApprox.filt;
Model = sysID.par.OL.time.initSys.prd.filt;
% Model = sysID.par.OL.time.fixedOrder1.prd.filt;
% Model = sysID.par.OL.time.fixedOrder5.prd.filt;
% Model = sysID.par.OL.time.straight.prd.filt;
% Model = G(1,1);

trainData   = sysID.CLdataFilt.train.direct;
trainData_U = sysID.CLdataFilt.train.u;
trainData_Y = sysID.CLdataFilt.train.y;

validData   = sysID.CLdataFilt.valid.direct;
validData_U = sysID.CLdataFilt.valid.u;
validData_Y = sysID.CLdataFilt.valid.y;
validData_T = sysID.CLdataFilt.valid.tVec;
% validData   = sysID.OLdataFilt.valid.id;
% validData_U = sysID.OLdataFilt.valid.u;
% validData_Y = sysID.OLdataFilt.valid.y;
% validData_T = sysID.OLdata.valid.tVec;

[test,~] = compare(trainData,Model,1);
[test_val,~] = compare(validData,Model,1);

predResiduals = trainData_Y-test.y;
predResiduals_val = validData_Y-test_val.y;


saveDataToNN(predResiduals,trainData_U)

saveFileNN = '_IDdata\NN\lstm_model_parallelData.keras';
epochs = 150;
rmseTrain = trainNN(saveFileNN,epochs);

[predNN_parll,errNN_parll] = predNN(saveFileNN,predResiduals_val,validData_U);

%%
figure(901);clf
    subplot(211)
        plot(predNN_parll);hold on
        plot(predResiduals_val);hold on
    subplot(212)
        plot(errNN_parll);hold on
        plot(predResiduals_val)
            legend([num2str(rms(errNN_parll))],[num2str(rms(predResiduals_val))])

Pred_trad = test_val.y;
Pred_parl = test_val.y+predNN_parll';

figure(902);clf
    subplot(211);hold on;grid minor
        plot(validData_T,Pred_trad)
        plot(validData_T,Pred_parl);
        plot(validData_T,validData_Y,'g--');
            xlabel('Time [s]')
            ylabel('Temperature [degC]')
            title('1-Step ahead temperature prediction of identified model')
            % legend('Using first order approx','First order + NN','Validation data','location','best')
            legend('Using initial system','Initial system + NN','Validation data','location','best')
            % legend('Using first principle model','first principle model + NN','Validation data','location','best')
            % legend('Using 5th order model','5th order model + NN','Validation data','location','best')
    subplot(212);hold on;grid minor
        plot(validData_T,Pred_trad-validData_Y)
        plot(validData_T,Pred_parl-validData_Y);
        yline(0,'k--')
            xlabel('Time [s]')
            ylabel('Residuals [degC]')
            title('1-Step ahead temperature prediction of identified model')
            legend([num2str(rms(Pred_trad-validData_Y))],[num2str(rms(Pred_parl-validData_Y))])
            legend(['Parametric, RMSE: ',num2str(rms(Pred_trad-validData_Y))],['Combination, RMSE: ',num2str(rms(Pred_parl-validData_Y))],'location','best')

%% Neural network training on residuals using closed loop models

% Model = sysID.par.CL.time.direct.filt;
% Model = sysID.par.CL.time.classic.filt;
% Model = sysID.par.CL.time.coprime.filt;
Model = sysID.par.CL.time.twostage.filt;
% Model = G(1,1);

% trainData   = sysID.CLdataFilt.train.direct;
% trainData_U = sysID.CLdataFilt.train.u;
% trainData_Y = sysID.CLdataFilt.train.y;
trainData   = sysID.OLdataFilt.train.id;
trainData_U = sysID.OLdataFilt.train.u;
trainData_Y = sysID.OLdataFilt.train.y;

% validData   = sysID.CLdataFilt.valid.direct;
% validData_U = sysID.CLdataFilt.valid.u;
% validData_Y = sysID.CLdataFilt.valid.y;
% validData_T = sysID.CLdataFilt.valid.tVec;
validData   = sysID.OLdataFilt.valid.id;
validData_U = sysID.OLdataFilt.valid.u;
validData_Y = sysID.OLdataFilt.valid.y;
validData_T = sysID.OLdata.valid.tVec;

[test,~] = compare(trainData,Model,1);
[test_val,~] = compare(validData,Model,1);

predResiduals = trainData_Y-test.y;
predResiduals_val = validData_Y-test_val.y;


saveDataToNN(predResiduals,trainData_U)

saveFileNN = '_IDdata\NN\lstm_model_parallelData.keras';
epochs = 150;
rmseTrain = trainNN(saveFileNN,epochs);

[predNN_parll,errNN_parll] = predNN(saveFileNN,predResiduals_val,validData_U);

%%
% figure(901);clf
%     subplot(211)
%         plot(predNN_parll);hold on
%         plot(predResiduals_val);hold on
%     subplot(212)
%         plot(errNN_parll);hold on
%         plot(predResiduals_val)
%             legend([num2str(rms(errNN_parll))],[num2str(rms(predResiduals_val))])

Pred_trad = test_val.y;
Pred_parl = test_val.y+predNN_parll';

figure(902);clf
    subplot(211);hold on;grid minor
        plot(validData_T,Pred_trad)
        plot(validData_T,Pred_parl);
        plot(validData_T,validData_Y,'g--');
            xlabel('Time [s]')
            ylabel('Temperature [degC]')
            title('1-Step ahead temperature prediction of identified model')
            % legend('Using CL, direct',' CL, direct + NN','Validation data','location','best')
            legend('Using two-stage model','Two-stage model + NN','Validation data','location','best')
            % legend('Using 5th order model','5th order model + NN','Validation data','location','best')
    subplot(212);hold on;grid minor
        plot(validData_T,Pred_trad-validData_Y)
        plot(validData_T,Pred_parl-validData_Y);
        yline(0,'k--')
            xlabel('Time [s]')
            ylabel('Residuals [degC]')
            title('1-Step ahead temperature prediction of identified model')
            legend([num2str(rms(Pred_trad-validData_Y))],[num2str(rms(Pred_parl-validData_Y))])
            legend(['Parametric, RMSE: ',num2str(rms(Pred_trad-validData_Y))],['Combination, RMSE: ',num2str(rms(Pred_parl-validData_Y))],'location','best')



































































