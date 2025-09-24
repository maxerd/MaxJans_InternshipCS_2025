clear all
% close all
clc

%% Main variables
baseFig_sim = 8000;

makeValSet  = true;
removeTrans = false;

%% Make the system used for data generation
% addpath(genpath('../..'))
disp('Constructing system model')
main_lumpedSystem_water_V3

G = sysTMC(32,[3 6]);

opt = bodeoptions('cstprefs');
opt.PhaseWrapping = 'on';
opt.FreqUnits = 'Hz';

%% Define some simulation variables
fs = 1;     % [Hz]    Sampling frequency
Ts = 1/fs;  % [s]     Sampling time
% tMeas = 20; % [hours] Measurement time
tMeas = 8; % [hours] Measurement time

T0 = 23; % [degC] Inital Temperature
Tamb = 23; % [degC] Ambient Temperature

T_ref = 45; % [degC] Reference Temperature

N = tMeas*3600*fs; % [-] Measurement samples
tVec = linspace(fs,tMeas*3600,N); % [s] Time vector

%% Define some controller variables
controllerBW = 0.001; % [Hz] Controller bandwidth

%% Define the noise on the signal
v     = 0.15.*randn(N,1);
v_val = 0.15.*randn(N,1);
% v = 0.*randn(N,1);

%% Define the input signal
maxAmplitude = 24; % [W] Maximum exitation signal
positiveOnly = 1;  % [-] Define whether the exitation signal can be negative

disp('Generating multisine for identification')
% dist = genMultisine(fs, N, 0.1, maxAmplitude, positiveOnly);
dist     = genMultisine(fs, N, 1, maxAmplitude,     positiveOnly);
dist_val = dist;%genMultisine(fs, N, 1, maxAmplitude/1.5, positiveOnly);
% dist = maxAmplitude.*ones(1,N);

%% Adapt the system for use in closed loop simulation
% Define the plant's in and outputnames
    G.inputName = {'HT3','Tamb'};
    G.outputName = 'TM1';

% Generate a controller
    K = genPID(G(1,1),controllerBW);
% Define the plant's in and outputnames
    K.tot.inputName = 'TM1_err';
    K.tot.outputName = 'cOut';

% Sums for interconnection of CL system
    refSum = sumblk('TM1_err = TM1_ref - TM1');
    distSum = sumblk('HT3 = cOut + dist');

% Define the total CL plant model
    P = connect(G,K.tot,refSum,distSum,{'TM1_ref','dist','Tamb'},{'TM1','HT3'});

%% Plot the noise signal
% Make the PSD of the noise and disturbance signals
    [amp_V,cas_V,freq_V] = fftCas_V2(v,fs);
    [amp_D,cas_D,freq_D] = fftCas_V2(dist,fs);

% Define FRD models of the PSD for easy use/plotting
    V      = frd(amp_V,freq_V);
    D      = frd(amp_D,freq_D);

% Visualize
    figure(baseFig_sim);clf
        subplot(211);
            semilogx(freq_V,db(amp_V));grid minor
                xlabel('Frequency [Hz]')
                ylabel('Power [degC^2/Hz]')
                title('Power of noise signal applied to simulation')
        subplot(212)
            semilogx(freq_D,db(amp_D));grid minor
                xlabel('Frequency [Hz]')
                ylabel('Power [degC^2/Hz]')
                title('Power of noise signal applied to simulation')
    
    disp(['Average power of noise signal: ',num2str(db(mean(amp_V))),'dB'])
    disp(['Maximal theoretically achievable accuracy (non-filtered data): ', num2str((1-mean(amp_V))*100),'%  ?????????????'])
    
    figure(baseFig_sim+1);
            bode(K.tot,opt);grid minor;

%% Simulate the systems

disp('Running simulation')

% Ambient temperature vector, constant ambient temperature for simulations
    ambVec = Tamb.*ones(1,length(dist));

% Define the initial temperatures (states) for the simulation, also constant.
    x0_OL = T0.*ones(size(G.A,1),1);

% Run the open loop simulation, with the previously defined noise signal added
    y_OL     = lsim(G,[dist    ;ambVec],tVec,x0_OL)+v;
    y_OL_val = lsim(G,[dist_val;ambVec],tVec,x0_OL)+v_val;


%% Plot the simulation results

lw = 1.5; % Define the wanted linewidth for the plots

figure(baseFig_sim+101);clf
    plot(tVec,y_OL,LineWidth=lw);grid minor
        xlabel('Time [s]')
        ylabel('Temperature [degC]')
        title(['Thermal mass temperature, open loop identification, ',num2str(tMeas),' hours'])

% figure(baseFig_sim+102);clf;hold on
%     yline(45,'r--',LineWidth=lw)
%     plot(tVec,y_CL(:,1),LineWidth=lw);grid minor
%         xlabel('Time [s]')
%         ylabel('Temperature [degC]')
%         title(['Thermal mass temperature, closed loop identification, ',num2str(tMeas),' hours'])
%         legend('Setpoint','Thermal mass temperature',Location='best')

%% Define the in and output vector sizes

% Define the range of the data that is used for identification and validation
if removeTrans
    idx1_OL = 8*3600*fs;
    idx1_CL = 1*3600*fs;
else
    idx1_OL = 1;
    idx1_CL = 1;
end
if makeValSet
    idx2_OL = N-(N-8*3600*fs)*0.2;
    idx2_CL = N-(N-1*3600*fs)*0.2;
else
    idx2_OL = N;
    idx2_CL = N;
end

idxRange_OL_trns = 1:idx1_OL+1;
idxRange_OL      = idx1_OL:idx2_OL;
idxRange_OL_val  = (idx2_OL-1):N;

% idxRange_CL_trns = 1:idx1_CL+1;
% idxRange_CL      = idx1_CL:idx2_CL;
% idxRange_CL_val  = (idx2_CL-1):N;

%% Define the in and output vectors, for easy and consistent use
% Open loop
outputData_OL_trns = y_OL(idxRange_OL_trns)-ambVec(idxRange_OL_trns)';
inputData_OL_trns  = dist(idxRange_OL_trns);

outputData_OL = y_OL(idxRange_OL)-ambVec(idxRange_OL)';
inputData_OL  = dist(idxRange_OL);

outputData_OL_val = y_OL_val-ambVec';
inputData_OL_val  = dist_val;

outputData_OL_full = y_OL-ambVec';
inputData_OL_full  = dist;

% % Closed loop for traditional non-parametric
% outputData_CL_trd = y_CL(idxRange_CL,2);
% inputData_CL_trd  = dist(idxRange_CL);
% 
% % Closed loop for LPM non-parametric
% outputData_CL_LPM = y_CL(idxRange_CL,1)-ambVec(idxRange_CL)';
% inputData_CL_LPM  = dist(idxRange_CL);
% 
% outputData_CL_full = y_CL(:,1)-ambVec';
% inputData_CL_full  = dist;

% Define the time vectors
tVec_OL_trns = tVec(idxRange_OL_trns);
tVec_OL      = tVec(idxRange_OL);
tVec_OL_val  = tVec;
tVec_OL_full = tVec;

% tVec_CL      = tVec(idxRange_CL);

%% Filter the data vectors to remove some of the effect of noise
% High order lowpass filter, to filter away all the noise above the 
% multisine exitation frequency
    LPfilt = orderLP(0.01,2,0,0,0.7)^2;
    % [LPfiltB,LPfiltA] = butter(4,0.005);
    % LPfilt = tf(LPfiltB,LPfiltA,1);

% Actually filter the data
    % Open loop data
    outputData_OL_trns_filt = lsim(LPfilt,outputData_OL_trns-outputData_OL_trns(1),tVec_OL_trns)'+outputData_OL_trns(1);
    inputData_OL_trns_filt = lsim(LPfilt,inputData_OL_trns-inputData_OL_trns(1),tVec_OL_trns)'+inputData_OL_trns(1);

    outputData_OL_filt = lsim(LPfilt,outputData_OL-outputData_OL(1),tVec_OL)'+outputData_OL(1);
    inputData_OL_filt = lsim(LPfilt,inputData_OL-inputData_OL(1),tVec_OL)'+inputData_OL(1);

    outputData_OL_val_filt = lsim(LPfilt,outputData_OL_val-outputData_OL_val(1),tVec_OL_val)'+outputData_OL_val(1);
    inputData_OL_val_filt = lsim(LPfilt,inputData_OL_val-inputData_OL_val(1),tVec_OL_val)'+inputData_OL_val(1);

    outputData_OL_full_filt = lsim(LPfilt,outputData_OL_full-outputData_OL_full(1),tVec_OL_full)'+outputData_OL_full(1);
    inputData_OL_full_filt = lsim(LPfilt,inputData_OL_full-inputData_OL_full(1),tVec_OL_full)'+inputData_OL_full(1);
    
    % % Closed loop data
    % outputData_CL_trd_filt = lsim(LPfilt,outputData_CL_trd-outputData_CL_trd(1),tVec_CL)'+outputData_CL_trd(1);
    % inputData_CL_trd_filt = lsim(LPfilt,inputData_CL_trd-inputData_CL_trd(1),tVec_CL)'+inputData_CL_trd(1);
    % 
    % outputData_CL_LPM_filt = lsim(LPfilt,outputData_CL_LPM-outputData_CL_LPM(1),tVec_CL)'+outputData_CL_LPM(1);
    % inputData_CL_LPM_filt = lsim(LPfilt,inputData_CL_LPM-inputData_CL_LPM(1),tVec_CL)'+inputData_CL_LPM(1);

%% FRF measurements using the new state vector
% Make the FRF's using the raw data
    [sysID.nonPar.trd.OL.raw, ~] = makeOpenLoopFRF_sysIdent(outputData_OL, inputData_OL(1,:), fs);
    % [sysID.nonPar.trd.CL.raw, ~] = makeClosedLoopFRF(outputData_CL_trd, inputData_CL_trd(1,:), K.tot, fs);

% Make the FRF's using the filtered data
    [sysID.nonPar.trd.OL.filt, ~] = makeOpenLoopFRF_sysIdent(outputData_OL_filt, inputData_OL_filt(1,:), fs);
    % [sysID.nonPar.trd.CL.filt, ~] = makeClosedLoopFRF(outputData_CL_trd_filt, inputData_CL_trd_filt(1,:), K.tot, fs);

% Define and do the non-parametric identification using the LPM, using both
% OL and CL data
    polyOrder = 6; % The order of the polynomial that is fitted to the data
    locality  = 8; % Amount of points (pos & neg) to consider around apprx freq
    
        disp('Starting OL LPM identification')
    [sysID.nonPar.lpm.OL.raw,~]  = sysID_LPM(outputData_OL     ,inputData_OL(1,:)     ,fs,locality,polyOrder);
    [sysID.nonPar.lpm.OL.filt,~] = sysID_LPM(outputData_OL_filt,inputData_OL_filt(1,:),fs,locality,polyOrder);
        disp('Finished OL LPM identification')
    %     disp('Starting CL LPM identification')
    % [sysID.nonPar.lpm.CL.raw,~]  = sysID_LPM(outputData_CL_LPM     ,inputData_CL_LPM(1,:)     ,fs,locality,polyOrder);
    % [sysID.nonPar.lpm.CL.filt,~] = sysID_LPM(outputData_CL_LPM_filt,inputData_CL_LPM_filt(1,:),fs,locality,polyOrder);
    %     disp('Finished CL LPM identification')

%% Visualization of the non-parametric identification
% Open loop data
    figure(baseFig_sim+201);clf
        subplot(211)
            bode(sysID.nonPar.trd.OL.raw,sysID.nonPar.lpm.OL.raw,G(1,1),'g--',opt);grid minor
            xlim([1e-5 1])
                xlabel('Frequency [Hz]')
                ylabel('Magnitude [dB]')
                title(['Heater to thermal mass temperature, open loop identification using raw data, ',num2str(tMeas),' hours'])
                legend('Traditional','LPM','Model')
        subplot(212)
            bode(sysID.nonPar.trd.OL.filt,sysID.nonPar.lpm.OL.filt,G(1,1),'g--',opt);grid minor
            xlim([1e-5 1])
                xlabel('Frequency [Hz]')
                ylabel('Magnitude [dB]')
                title(['Heater to thermal mass temperature, open loop identification using filtered data, ',num2str(tMeas),' hours'])
                legend('Traditional','LPM','Model')

% % Closed loop data
%     figure(baseFig_sim+202);clf
%         subplot(211);
%             bode(sysID.nonPar.trd.CL.raw,sysID.nonPar.lpm.CL.raw,G(1,1),'g--',opt);hold on;grid minor
%             xline(controllerBW,'m--')
%             xlim([1e-5 1])
%                 xlabel('Frequency [Hz]')
%                 ylabel('Magnitude [dB]')
%                 title(['Heater to thermal mass temperature, closed loop identification using raw data, ',num2str(tMeas),' hours'])
%                 legend('Traditional','LPM','Model','Controller bandwidth')
%         subplot(212);
%             bode(sysID.nonPar.trd.CL.filt,sysID.nonPar.lpm.CL.filt,G(1,1),'g--',opt);hold on;grid minor
%             xline(controllerBW,'m--')
%             xlim([1e-5 1])
%                 xlabel('Frequency [Hz]')
%                 ylabel('Magnitude [dB]')
%                 title(['Heater to thermal mass temperature, closed loop identification using filtered data, ',num2str(tMeas),' hours'])
%                 legend('Traditional','LPM','Model','Controller bandwidth')

%% First order parametric approximation, using time data
    sysID.par.firstAprrox.OL.raw  = step_sysID(inputData_OL_full',zeros(size(inputData_OL_full))',outputData_OL_full,tVec_OL_full,10000*Ts);
    sysID.par.firstAprrox.OL.filt = step_sysID(inputData_OL_full_filt',zeros(size(inputData_OL_full_filt))',outputData_OL_full_filt',tVec_OL_full,10000*Ts);

%% Visualization of the first order approximation
figure(baseFig_sim+203);clf
    bode(sysID.par.firstAprrox.OL.raw,sysID.par.firstAprrox.OL.filt,G(1,1),'g--',opt);grid minor
    xlim([1e-5 0.1])
        xlabel('Frequency [Hz]')
        ylabel('Magnitude [dB]')
        title(['Heater to thermal mass temperature, first order approximation, ',num2str(tMeas),' hours'])
        legend('Unfiltered data','Filtered data','Model')

%% Parametric system identification, using time data, pre-requisites

% Definitions for identification
    nx = 4; % Model order for fixed order identification

% nx order initial system
    % init_sys = idss([sysID.par.firstAprrox.OL.raw.A 0 0;0 sysID.par.firstAprrox.OL.raw.A*10 0;0 0 sysID.par.firstAprrox.OL.raw.A*100],[sysID.par.firstAprrox.OL.raw.B;0;0],[sysID.par.firstAprrox.OL.raw.C 0 0],0,zeros(3,1),zeros(3,1),0);
    initA = sysID.par.firstAprrox.OL.raw.A.*eye(nx);
    initB = zeros(nx,1);initB(1) = sysID.par.firstAprrox.OL.raw.B;
    initC = zeros(1,nx);initC(1) = sysID.par.firstAprrox.OL.raw.C;
    
    init_sys = idss(initA,initB,initC,0,zeros(nx,1),zeros(nx,1),0);
    init_sys.Structure.B.Free = ones(nx,1);
    init_sys.Structure.C.Free = ones(1,nx);
    init_sys.Structure.K.Free = zeros(nx,1);

% Define the data as iddata's
    OL_dat      = iddata(outputData_OL          ,inputData_OL(:,:)'    ,Ts);
    OL_dat_filt = iddata(outputData_OL_filt'    ,inputData_OL_filt'    ,Ts);
    % CL_dat      = iddata(outputData_CL_LPM      ,inputData_CL_LPM(:,:)',Ts);
    % CL_dat_filt = iddata(outputData_CL_LPM_filt',inputData_CL_LPM_filt',Ts);

% Define the data that is not used in the identification as validation data

if makeValSet
    OL_dat_val      = iddata(outputData_OL_val      ,inputData_OL_val',Ts);
    OL_dat_val_filt = iddata(outputData_OL_val_filt',inputData_OL_val',Ts);
else
    OL_dat_val      = iddata(1,1,Ts);
    OL_dat_val_filt = iddata(1,1,Ts);
end
if removeTrans
    OL_dat_trns      = iddata(outputData_OL_trns      ,inputData_OL_trns',Ts);
    OL_dat_trns_filt = iddata(outputData_OL_trns_filt',inputData_OL_trns',Ts);
else
    OL_dat_trns      = iddata(1,1,Ts);
    OL_dat_trns_filt = iddata(1,1,Ts);
end

optSS_sim = ssestOptions('Focus','Simulation');
optSS_prd = ssestOptions('Focus','Prediction');

%% Parametric system identification, using time data

% % All tested parametric identification options using OL data, using simulation focus
%         disp('Open Loop identification using an initial system: in progress')
%     sysID.par.initSys.sim.OL.raw  = ssest(OL_dat,init_sys,optSS_sim);
%     sysID.par.initSys.sim.OL.filt = ssest(OL_dat_filt,init_sys,optSS_sim);
%         disp('Open Loop identification using an initial system: done')
% 
%         disp(['Open Loop identification using a fixed order (nx=',num2str(nx),'): in progress'])
%     sysID.par.fixedOrder.sim.OL.raw  = ssest(OL_dat,nx,optSS_sim);
%     sysID.par.fixedOrder.sim.OL.filt = ssest(OL_dat_filt,nx,optSS_sim);
%         disp(['Open Loop identification using a fixed order (nx=',num2str(nx),'): done'])
% 
%         disp('Open Loop identification using the raw data: in progress')
%     sysID.par.straight.sim.OL.raw  = ssest(OL_dat);
%     sysID.par.straight.sim.OL.filt = ssest(OL_dat_filt);
%         disp('Open Loop identification using the raw data: done')
% 
% % All tested parametric identification options using CL data, using simulation focus
%         disp('Closed Loop identification using an initial system: in progress')
%     sysID.par.initSys.sim.CL.raw  = ssest(CL_dat,init_sys,optSS_sim);
%     sysID.par.initSys.sim.CL.filt = ssest(CL_dat_filt,init_sys,optSS_sim);
%         disp('Closed Loop identification using an initial system: done')
% 
%         disp(['Closed Loop identification using a fixed order (nx=',num2str(nx),'): in progress'])
%     sysID.par.fixedOrder.sim.CL.raw  = ssest(CL_dat,nx,optSS_sim);
%     sysID.par.fixedOrder.sim.CL.filt = ssest(CL_dat_filt,nx,optSS_sim);
%         disp(['Closed Loop identification using a fixed order (nx=',num2str(nx),'): done'])
% 
%         disp('Closed Loop identification using filtered data: in progress')
%     sysID.par.straight.sim.CL.raw  = ssest(CL_dat);
%     sysID.par.straight.sim.CL.filt = ssest(CL_dat_filt);
%         disp('Closed Loop identification using filtered data: done')

% All tested parametric identification options using OL data, using prediction focus
        disp('Open Loop identification using an initial system: in progress')
    sysID.par.initSys.sim.OL.raw  = ssest(OL_dat,init_sys,optSS_prd);
    sysID.par.initSys.sim.OL.filt = ssest(OL_dat_filt,init_sys,optSS_prd);
        disp('Open Loop identification using an initial system: done')

        disp(['Open Loop identification using a fixed order (nx=',num2str(nx),'): in progress'])
    sysID.par.fixedOrder.sim.OL.raw  = ssest(OL_dat,nx,optSS_prd);
    sysID.par.fixedOrder.sim.OL.filt = ssest(OL_dat_filt,nx,optSS_prd);
        disp(['Open Loop identification using a fixed order (nx=',num2str(nx),'): done'])

        disp('Open Loop identification using the raw data: in progress')
    sysID.par.straight.sim.OL.raw  = ssest(OL_dat);
    sysID.par.straight.sim.OL.filt = ssest(OL_dat_filt);
        disp('Open Loop identification using the raw data: done')

% % All tested parametric identification options using CL data, using simulation focus
%         disp('Closed Loop identification using an initial system: in progress')
%     sysID.par.initSys.sim.CL.raw  = ssest(CL_dat,init_sys,optSS_prd);
%     sysID.par.initSys.sim.CL.filt = ssest(CL_dat_filt,init_sys,optSS_prd);
%         disp('Closed Loop identification using an initial system: done')
% 
%         disp(['Closed Loop identification using a fixed order (nx=',num2str(nx),'): in progress'])
%     sysID.par.fixedOrder.sim.CL.raw  = ssest(CL_dat,nx,optSS_prd);
%     sysID.par.fixedOrder.sim.CL.filt = ssest(CL_dat_filt,nx,optSS_prd);
%         disp(['Closed Loop identification using a fixed order (nx=',num2str(nx),'): done'])
% 
%         disp('Closed Loop identification using filtered data: in progress')
%     sysID.par.straight.sim.CL.raw  = ssest(CL_dat);
%     sysID.par.straight.sim.CL.filt = ssest(CL_dat_filt);
%         disp('Closed Loop identification using filtered data: done')

%% Visualize and compare the identifications in frequency domain

figure(baseFig_sim+301);clf
    bode(sysID.par.initSys.sim.OL.raw,opt);hold on
    bode(sysID.par.fixedOrder.sim.OL.raw,opt);hold on
    bode(sysID.par.firstAprrox.OL.raw,opt);hold on
    bode(sysID.par.straight.sim.OL.raw,opt);
    bode(sysTMC(32,3),'g--');grid minor
        title('Bode plots of identification using raw Open Loop data')
        legend('InitSysIdent','Fixed order','First order approx','Straight data','Model')

figure(baseFig_sim+302);clf
    bode(sysID.par.initSys.sim.OL.filt,opt);hold on
    bode(sysID.par.fixedOrder.sim.OL.filt,opt);hold on
    bode(sysID.par.firstAprrox.OL.filt,opt);hold on
    bode(sysID.par.straight.sim.OL.filt,opt);
    bode(sysTMC(32,3),'g--');grid minor
        title('Bode plots of identification using filtered Open Loop data')
        legend('InitSysIdent','Fixed order','First order approx','Straight data','Model')

% figure(baseFig_sim+303);clf
%     bode(sysID.par.initSys.sim.CL.raw,opt);hold on
%     bode(sysID.par.fixedOrder.sim.CL.raw,opt)
%     bode(sysID.par.straight.sim.CL.raw,opt);
%     bode(sysTMC(32,3),'g--');grid minor
%         title('Bode plots of identification using raw Closed Loop data')
%         legend('InitSysIdent','Fixed order','Straight data','Model')
% 
% figure(baseFig_sim+304);clf
%     bode(sysID.par.initSys.sim.CL.filt,opt);hold on
%     bode(sysID.par.fixedOrder.sim.CL.filt,opt)
%     bode(sysID.par.straight.sim.CL.filt,opt);
%     bode(sysTMC(32,3),'g--');grid minor
%         title('Bode plots of identification using filtered Closed Loop data')
%         legend('InitSysIdent','Fixed order','Straight data','Model')


figure(baseFig_sim+305);clf
    subplot(221)
        bode(sysID.par.initSys.sim.OL.raw,opt);hold on
        bode(sysID.par.initSys.sim.OL.filt,opt);hold on
        bode(sysTMC(32,3),'g--');grid minor
            title('Initial system')
            legend('Raw','Filtered')
    subplot(222)
        bode(sysID.par.fixedOrder.sim.OL.raw,opt);hold on
        bode(sysID.par.fixedOrder.sim.OL.filt,opt);hold on
        bode(sysTMC(32,3),'g--');grid minor
            title(['Fixed order, nx=',num2str(nx)])
            legend('Raw','Filtered')
    subplot(223)
        bode(sysID.par.firstAprrox.OL.raw,opt);hold on
        bode(sysID.par.firstAprrox.OL.filt,opt);hold on
        bode(sysTMC(32,3),'g--');grid minor
            title('First order approximation')
            legend('Raw','Filtered')
    subplot(224)
        bode(sysID.par.straight.sim.OL.raw,opt);hold on
        bode(sysID.par.straight.sim.OL.filt,opt);hold on
        bode(sysTMC(32,3),'g--');grid minor
            title('Straight data')
            legend('Raw','Filtered')

%% Validate and compare the identifications in time domain, using x-step ahead prediction and a residual test

% Choose which dataset to use for the validation
    % dataSet = OL_dat;
    dataSet = OL_dat_val;
    % dataSet = OL_dat_trns;
    
    % dataSet_filt = OL_dat_filt;
    dataSet_filt = OL_dat_val_filt;
    % dataSet_filt = OL_dat_trns_filt;

% Choose the amount of steps to use for the prediction (inf=simulation)
    % xStep = 1;
    xStep = inf;

% Run the predictions and plot them
    sysID_predTests_V3

% Run the residual tests and plot them
    sysID_resTests_V3


%% Frequency domain identification

nx_freq = 4;

OL_freqDat_trd      = idfrd(squeeze(sysID.nonPar.trd.OL.raw.ResponseData) ,sysID.nonPar.trd.OL.raw.Frequency ,1);
OL_freqDat_lpm      = idfrd(squeeze(sysID.nonPar.lpm.OL.raw.ResponseData) ,sysID.nonPar.lpm.OL.raw.Frequency ,1);
OL_freqDat_trd_filt = idfrd(squeeze(sysID.nonPar.trd.OL.filt.ResponseData),sysID.nonPar.trd.OL.filt.Frequency,1);
OL_freqDat_lpm_filt = idfrd(squeeze(sysID.nonPar.lpm.OL.filt.ResponseData),sysID.nonPar.lpm.OL.filt.Frequency,1);

sysID.par.freqTRD.sim.OL.raw = ssest(OL_freqDat_trd,nx_freq,optSS_sim);
sysID.par.freqLPM.sim.OL.raw = ssest(OL_freqDat_lpm,nx_freq,optSS_sim);
% sysID.par.freqTRD.sim.OL.raw  = ssest(OL_freqDat_trd     ,init_sys,optSS_sim);
% sysID.par.freqLPM.sim.OL.raw  = ssest(OL_freqDat_lpm     ,init_sys,optSS_sim);
% sysID.par.freqTRD.sim.OL.raw  = ssest(OL_freqDat_trd     ,sysID.par.fixedOrder.sim.OL.raw,optSS_sim);
% sysID.par.freqLPM.sim.OL.raw  = ssest(OL_freqDat_lpm     ,sysID.par.fixedOrder.sim.OL.raw,optSS_sim);
% sysID.par.freqTRD.sim.OL.filt = ssest(OL_freqDat_trd_filt,sysID.par.fixedOrder.sim.OL.raw,optSS_sim);
% sysID.par.freqLPM.sim.OL.filt = ssest(OL_freqDat_lpm_filt,sysID.par.fixedOrder.sim.OL.raw,optSS_sim);

%%
figure(baseFig_sim+701);clf
    % bode(sysID.par.freqTRD.sim.OL.raw,sysID.par.freqLPM.sim.OL.raw,sysID.par.freqTRD.sim.OL.filt,sysID.par.freqLPM.sim.OL.filt,G(1,1),'g--',opt);grid minor
    bode(sysID.par.freqTRD.sim.OL.raw,sysID.par.freqLPM.sim.OL.raw,G(1,1),'g--',opt);grid minor

% figure(baseFig_sim+801);clf
%     resid(OL_freqDat_trd,sysID.par.freqTRD.sim.OL.raw,sysID.par.freqLPM.sim.OL.raw);grid minor
% figure(baseFig_sim+802);clf
%     resid(OL_freqDat_lpm,sysID.par.freqTRD.sim.OL.raw,sysID.par.freqLPM.sim.OL.raw);grid minor

% figure(baseFig_sim+801);clf
%     resid(OL_dat,sysID.par.freqTRD.sim.OL.raw,sysID.par.freqLPM.sim.OL.raw);grid minor
% figure(baseFig_sim+802);clf
%     resid(OL_dat,sysID.par.freqTRD.sim.OL.raw,sysID.par.freqLPM.sim.OL.raw);grid minor

%%
dat.sysID.par.freqTRD.sim.OL.raw = lsim(sysID.par.freqTRD.sim.OL.raw,OL_dat.u,tVec_OL);

figure(baseFig_sim+901);clf
    plot(tVec_OL,dat.sysID.par.freqTRD.sim.OL.raw,tVec_OL,OL_dat.y)






























