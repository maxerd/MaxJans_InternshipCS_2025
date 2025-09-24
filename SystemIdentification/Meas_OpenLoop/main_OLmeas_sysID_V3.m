clear all
% close all
clc

%% Pre-requisites
opt = bodeoptions('cstprefs');
opt.PhaseWrapping = 'on';
opt.FreqUnits = 'Hz';

lw = 1.5; % Define the wanted linewidth for the plots

list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end

%% Main variables

% Figure number of of which all others are based
baseFig_sim = 8000;

% Some user inputs
makeValSet  = true;  % [true/false] Is a seperate validation set needed?
removeTrans = false; % [true/false] Does the transient need to be removed

%% Define the measurement paths
dataDir_OL = 'C:\Users\maxja\Documents\(4)School\Master\Q9_Internship\matlabFiles\measurements\processedData\p__CB_none__ST_ms__SM_24w__SA_12w__DT_250829__MD_8h__WT_no__DS_1000.mat'; % Open loop data used for identification
dataDir_CL = 'C:\Users\maxja\Documents\(4)School\Master\Q9_Internship\matlabFiles\measurements\processedData\p__CB_none__ST_ms__SM_24w__SA_12w__DT_250905__MD_9h__WT_no__DS_100.mat'; % Closed loop data used for validation

% Load in the measurement data
disp('Loading in measurement data...')
    dat   = load(dataDir_OL);
    datCL = load(dataDir_CL);
disp('Loading done!')

%% Define some measurement variables for easier use
fs = 1/dat.Ts;     % [Hz]    Sampling frequency
Ts = dat.Ts;  % [s]     Sampling time
N = size(dat.tVec,2); % [-] Measurement samples

tMeas = N*Ts/3600; % [hours] Measurement time

tVec = dat.tVec; % [s] Time vector

%% Define the input signal
dist     = dat.Watt./5;

%% Plot the noise signal
% Make the PSD of the noise and disturbance signals
    [amp_D,cas_D,freq_D] = fftCas_V2(dist,fs);

% Define FRD model of the PSD for easy use/plotting
    D      = frd(amp_D,freq_D);

% Visualize the disturbance signal's power
    figure(baseFig_sim);clf
        semilogx(freq_D,db(amp_D));grid minor
            xlabel('Frequency [Hz]')
            ylabel('Power [$degC^2/Hz$]')
            title('Power of disturbance/identification signal applied to system')
    
%% Simulate the systems

% Ambient temperature vector
    ambVec = dat.tempAmb;

% Define the system output, for easier use later
    y_OL   = dat.tempTM;

%% Plot the measurement results
figure(baseFig_sim+101);clf
    set(gcf,'position',[700 100 700 300])
    plot(tVec,y_OL,LineWidth=lw);grid minor
        xlabel('Time [s]')
        ylabel('Temperature [degC]')
        title(['Thermal mass temperature, open loop identification, ',num2str(tMeas,3),' hours'])

%% Define the range of the data that is used for identification and validation
if removeTrans
    idx1_OL = 8*3600*fs;
else
    idx1_OL = 1;
end
if makeValSet
    idx2_OL = round(N-(N-4*3600*fs)*0.2);
else
    idx2_OL = N;
end

idxRange_OL_trns = 1:idx1_OL+1;     % Transient data range
idxRange_OL      = idx1_OL:idx2_OL; % Training data range
idxRange_OL_val  = (idx2_OL-1):N;   % Validation data range

%% Define the in and output vectors, for easy and consistent use
sysID.data.OL.trans.out = y_OL(idxRange_OL_trns)-ambVec(idxRange_OL_trns);
sysID.data.OL.trans.in  = dist(idxRange_OL_trns);

sysID.data.OL.train.out = y_OL(idxRange_OL)-ambVec(idxRange_OL);
sysID.data.OL.train.in  = dist(idxRange_OL);

sysID.data.OL.val.out = y_OL(idxRange_OL_val)-ambVec(idxRange_OL_val);
sysID.data.OL.val.in  = dist(idxRange_OL_val);

sysID.data.OL.full.out = y_OL-ambVec;
sysID.data.OL.full.in  = dist;

% Define the time vectors
sysID.data.OL.trans.tVec = tVec(idxRange_OL_trns);
sysID.data.OL.train.tVec = tVec(idxRange_OL);
sysID.data.OL.val.tVec   = tVec(idxRange_OL_val);
sysID.data.OL.full.tVec  = tVec;

sysID.data.CL.full.tVec  = datCL.tVec;

%% Filter the data vectors to remove some of the effect of noise
% High order lowpass filter, to filter away all the noise above the 
% multisine exitation frequency
    LPfilt = orderLP(0.01,2,0,0,0.7)^2;

% Actually filter the data
    % Open loop data
    sysID.dataFilt.OL.trans.out = lsim(LPfilt,sysID.data.OL.trans.out-sysID.data.OL.trans.out(1),sysID.data.OL.trans.tVec)'+sysID.data.OL.trans.out(1);
    sysID.dataFilt.OL.trans.in  = lsim(LPfilt,sysID.data.OL.trans.in-sysID.data.OL.trans.in(1)  ,sysID.data.OL.trans.tVec)'+sysID.data.OL.trans.in(1);

    sysID.dataFilt.OL.train.out = lsim(LPfilt,sysID.data.OL.train.out-sysID.data.OL.train.out(1),sysID.data.OL.train.tVec)'+sysID.data.OL.train.out(1);
    sysID.dataFilt.OL.train.in  = lsim(LPfilt,sysID.data.OL.train.in-sysID.data.OL.train.in(1)  ,sysID.data.OL.train.tVec)'+sysID.data.OL.train.in(1);

    sysID.dataFilt.OL.val.out   = lsim(LPfilt,sysID.data.OL.val.out-sysID.data.OL.val.out(1)    ,sysID.data.OL.val.tVec)'  +sysID.data.OL.val.out(1);
    sysID.dataFilt.OL.val.in    = lsim(LPfilt,sysID.data.OL.val.in-sysID.data.OL.val.in(1)      ,sysID.data.OL.val.tVec)'  +sysID.data.OL.val.in(1);

    sysID.dataFilt.OL.full.out  = lsim(LPfilt,sysID.data.OL.full.out-sysID.data.OL.full.out(1)  ,sysID.data.OL.full.tVec)' +sysID.data.OL.full.out(1);
    sysID.dataFilt.OL.full.in   = lsim(LPfilt,sysID.data.OL.full.in-sysID.data.OL.full.in(1)    ,sysID.data.OL.full.tVec)' +sysID.data.OL.full.in(1);
    
%% FRF measurements using the new state vector
% Make the FRF's using the raw data
    [sysID.nonPar.trd.OL.raw, ~] = makeOpenLoopFRF_sysIdent(sysID.data.OL.train.out, sysID.data.OL.train.in(1,:), fs);

% Make the FRF's using the filtered data
    [sysID.nonPar.trd.OL.filt, ~] = makeOpenLoopFRF_sysIdent(sysID.dataFilt.OL.train.out, sysID.dataFilt.OL.train.in(1,:), fs);

% Define and do the non-parametric identification using the LPM, using both
% OL and CL data
    polyOrder = 6; % The order of the polynomial that is fitted to the data
    locality  = 8; % Amount of points (pos & neg) to consider around apprx freq
    
        disp('Starting OL LPM identification')
    [sysID.nonPar.lpm.OL.raw,~]  = sysID_LPM(sysID.data.OL.train.out     ,sysID.data.OL.train.in(1,:)     ,fs,locality,polyOrder);
    [sysID.nonPar.lpm.OL.filt,~] = sysID_LPM(sysID.dataFilt.OL.train.out,sysID.dataFilt.OL.train.in(1,:),fs,locality,polyOrder);
        disp('Finished OL LPM identification')

%% Visualization of the non-parametric identification
% Open loop data
    figure(baseFig_sim+201);clf
    set(gcf,'position',[500 100 900 500])
    sgtitle(['Heater to TM temperature, open loop identification using unfiltered and filtered data, ',num2str(tMeas),' hours'])
        subplot(221)
            simpleBodemag(sysID.nonPar.trd.OL.raw,'Hz',lw);hold on;grid minor
            simpleBodemag(sysID.nonPar.lpm.OL.raw,'Hz',lw);grid minor
            xlim([1e-4 0.01])
                title(['Using unfiltered data'])
                legend('Traditional','LPM','Model','location','best')
        subplot(223)
            simpleBodephase(sysID.nonPar.trd.OL.raw,'Hz',lw,'wrap');hold on;grid minor
            simpleBodephase(sysID.nonPar.lpm.OL.raw,'Hz',lw,'wrap');grid minor
            xlim([1e-4 0.01])
                legend('Traditional','LPM','Model','location','best')
        subplot(222)
            simpleBodemag(sysID.nonPar.trd.OL.filt,'Hz',lw);hold on;grid minor
            simpleBodemag(sysID.nonPar.lpm.OL.filt,'Hz',lw);grid minor
            xlim([1e-4 0.01])
                title(['Using filtered data'])
                legend('Traditional','LPM','Model','location','best')
        subplot(224)
            simpleBodephase(sysID.nonPar.trd.OL.filt,'Hz',lw,'wrap');hold on;grid minor
            simpleBodephase(sysID.nonPar.lpm.OL.filt,'Hz',lw,'wrap');grid minor
            xlim([1e-4 0.01])
                legend('Traditional','LPM','Model','location','best')

%% First order parametric approximation, using time data
    sysID.par.firstAprrox.OL.raw  = step_sysID(sysID.data.OL.full.in',zeros(size(sysID.data.OL.full.in))',sysID.data.OL.full.out,sysID.data.OL.full.tVec,10000*fs);
    sysID.par.firstAprrox.OL.filt = step_sysID(sysID.dataFilt.OL.full.in',zeros(size(sysID.dataFilt.OL.full.in))',sysID.dataFilt.OL.full.out',sysID.data.OL.full.tVec,10000*fs);

%% Visualization of the first order approximation
figure(baseFig_sim+203);clf
    bode(sysID.par.firstAprrox.OL.raw,sysID.par.firstAprrox.OL.filt,opt);grid minor
    xlim([1e-5 0.1])
        xlabel('Frequency [Hz]')
        ylabel('Magnitude [dB]')
        title(['Heater to thermal mass temperature, first order approximation, ',num2str(tMeas),' hours'])
        legend('Unfiltered data','Filtered data','Model')

%% Parametric system identification, using time data, pre-requisites

% Definitions for identification
    nx = 4; % Model order for fixed order identification

% nx order initial system
    initA = sysID.par.firstAprrox.OL.raw.A.*eye(nx);
    initB = zeros(nx,1);initB(1) = sysID.par.firstAprrox.OL.raw.B;
    initC = zeros(1,nx);initC(1) = sysID.par.firstAprrox.OL.raw.C;
    
    init_sys = idss(initA,initB,initC,0,zeros(nx,1),zeros(nx,1),0);
    init_sys.Structure.B.Free = ones(nx,1);
    init_sys.Structure.C.Free = ones(1,nx);
    init_sys.Structure.K.Free = zeros(nx,1);
    
% Define the data as iddata's
    OL_dat      = iddata(sysID.data.OL.train.out'     ,sysID.data.OL.train.in(:,:)',Ts);
    CL_dat      = iddata((datCL.tempTM-datCL.tempAmb)',datCL.Watt'                 ,datCL.Ts);

    OL_dat_filt = iddata(sysID.dataFilt.OL.train.out' ,sysID.dataFilt.OL.train.in' ,Ts);

% Define the data that is not used in the identification as validation data
if makeValSet
    OL_dat_val      = iddata(sysID.data.OL.val.out'      ,sysID.data.OL.val.in',Ts);
    OL_dat_val_filt = iddata(sysID.dataFilt.OL.val.out',sysID.dataFilt.OL.val.in',Ts);
else
    OL_dat_val      = iddata(1,1,Ts);
    OL_dat_val_filt = iddata(1,1,Ts);
end
if removeTrans
    OL_dat_trns      = iddata(sysID.data.OL.trans.out'      ,sysID.data.OL.trans.in',Ts);
    OL_dat_trns_filt = iddata(sysID.dataFilt.OL.trans.out',sysID.dataFilt.OL.trans.in',Ts);
else
    OL_dat_trns      = iddata(1,1,Ts);
    OL_dat_trns_filt = iddata(1,1,Ts);
end

optSS_sim = ssestOptions('Focus','Simulation');
optSS_prd = ssestOptions('Focus','Prediction');

%% Parametric system identification, using time data

% All tested parametric identification options using OL data, using prediction focus
        disp('Open Loop identification using an initial system: in progress')
    sysID.par.initSys.sim.OL.raw  = ssest(OL_dat,init_sys,optSS_sim);
    sysID.par.initSys.sim.OL.filt = ssest(OL_dat_filt,init_sys,optSS_sim);
        disp('Open Loop identification using an initial system: done')

        disp(['Open Loop identification using a fixed order (nx=',num2str(nx),'): in progress'])
    sysID.par.fixedOrder.sim.OL.raw  = ssest(OL_dat,nx,optSS_sim);
    sysID.par.fixedOrder.sim.OL.filt = ssest(OL_dat_filt,nx,optSS_sim);
        disp(['Open Loop identification using a fixed order (nx=',num2str(nx),'): done'])

% % All tested parametric identification options using OL data, using prediction focus
%         disp('Open Loop identification using an initial system: in progress')
%     sysID.par.initSys.sim.OL.raw  = ssest(OL_dat,init_sys,optSS_prd);
%     sysID.par.initSys.sim.OL.filt = ssest(OL_dat_filt,init_sys,optSS_prd);
%         disp('Open Loop identification using an initial system: done')
% 
%         disp(['Open Loop identification using a fixed order (nx=',num2str(nx),'): in progress'])
%     sysID.par.fixedOrder.sim.OL.raw  = ssest(OL_dat,nx,optSS_prd);
%     sysID.par.fixedOrder.sim.OL.filt = ssest(OL_dat_filt,nx,optSS_prd);
%         disp(['Open Loop identification using a fixed order (nx=',num2str(nx),'): done'])


        disp('Open Loop identification using the raw data: in progress')
    sysID.par.straight.sim.OL.raw  = ssest(OL_dat);
    sysID.par.straight.sim.OL.filt = ssest(OL_dat_filt);
        disp('Open Loop identification using the raw data: done')

%% Visualize and compare the identifications in frequency domain

figure(baseFig_sim+301);clf
    bode(sysID.par.initSys.sim.OL.raw,opt);hold on
    bode(sysID.par.fixedOrder.sim.OL.raw,opt);hold on
    bode(sysID.par.firstAprrox.OL.raw,opt);hold on
    bode(sysID.par.straight.sim.OL.raw,opt);
        title('Bode plots of identification using raw Open Loop data')
        legend('InitSysIdent','Fixed order','First order approx','Straight data')

figure(baseFig_sim+302);clf
    bode(sysID.par.initSys.sim.OL.filt,opt);hold on
    bode(sysID.par.fixedOrder.sim.OL.filt,opt);hold on
    bode(sysID.par.firstAprrox.OL.filt,opt);hold on
    bode(sysID.par.straight.sim.OL.filt,opt);
        title('Bode plots of identification using filtered Open Loop data')
        legend('InitSysIdent','Fixed order','First order approx','Straight data')

figure(baseFig_sim+305);clf
    subplot(221)
        bode(sysID.par.initSys.sim.OL.raw,opt);hold on
        bode(sysID.par.initSys.sim.OL.filt,opt);hold on
            title('Initial system')
            legend('Raw','Filtered')
    subplot(222)
        bode(sysID.par.fixedOrder.sim.OL.raw,opt);hold on
        bode(sysID.par.fixedOrder.sim.OL.filt,opt);hold on
            title(['Fixed order, nx=',num2str(nx)])
            legend('Raw','Filtered')
    subplot(223)
        bode(sysID.par.firstAprrox.OL.raw,opt);hold on
        bode(sysID.par.firstAprrox.OL.filt,opt);hold on
            title('First order approximation')
            legend('Raw','Filtered')
    subplot(224)
        bode(sysID.par.straight.sim.OL.raw,opt);hold on
        bode(sysID.par.straight.sim.OL.filt,opt);hold on
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
    % xStep = 10;
    xStep = inf;

% Run the predictions and plot them
    OLmeas_sysID_predTests_V3

% Run the residual tests and plot them
    OLmeas_sysID_resTests_V3

%% Time based validation using different (CL) data set
    
err_initSys = lsim(sysID.par.initSys.sim.OL.filt,datCL.Watt,datCL.tVec)+datCL.tempAmb'-datCL.tempTM';
err_fixedOrder = lsim(sysID.par.fixedOrder.sim.OL.filt,datCL.Watt,datCL.tVec)+datCL.tempAmb'-datCL.tempTM';
err_firstAprrox = lsim(sysID.par.firstAprrox.OL.filt,datCL.Watt,datCL.tVec)+datCL.tempAmb'-datCL.tempTM';
err_straight = lsim(sysID.par.straight.sim.OL.filt,datCL.Watt,datCL.tVec)+datCL.tempAmb'-datCL.tempTM';
    
figure(baseFig_sim+651);clf
    subplot(211);hold on;grid minor
        plot(datCL.tVec,lsim(sysID.par.initSys.sim.OL.filt,datCL.Watt,datCL.tVec)+datCL.tempAmb')
        plot(datCL.tVec,lsim(sysID.par.fixedOrder.sim.OL.filt,datCL.Watt,datCL.tVec)+datCL.tempAmb')
        plot(datCL.tVec,lsim(sysID.par.firstAprrox.OL.filt,datCL.Watt,datCL.tVec)+datCL.tempAmb')
        plot(datCL.tVec,lsim(sysID.par.straight.sim.OL.filt,datCL.Watt,datCL.tVec)+datCL.tempAmb')
        plot(datCL.tVec,datCL.tempTM)
            xlabel('Time [s]')
            ylabel('Temperature [degC]')
            title('Comparison of measured and simulated data using identified model')
            legend('Simulated data using identified model (initSys)', ...
                   'Simulated data using identified model (fixed order)', ...
                   'Simulated data using identified model (firstApprox)', ...
                   'Simulated data using identified model (Straight)', ...
                   'Measured data','location','southeast')
    subplot(212);hold on;grid minor
        plot(datCL.tVec,err_initSys)
        plot(datCL.tVec,err_fixedOrder)
        plot(datCL.tVec,err_firstAprrox)
        plot(datCL.tVec,err_straight)
            xlabel('Time [s]')
            ylabel('Temperature diff. [degC]')
            title('Difference between measured and simulated data using identified model')
            legend(['Simulated data using identified model (initSys), RMSE: ',num2str(rms(err_initSys),2),'degC'], ...
                   ['Simulated data using identified model (fixed order), RMSE: ',num2str(rms(err_fixedOrder),2),'degC'], ...
                   ['Simulated data using identified model (firstApprox), RMSE: ',num2str(rms(err_firstAprrox),2),'degC'], ...
                   ['Simulated data using identified model (Straight), RMSE: ',num2str(rms(err_straight),2),'degC'], ...
               'location','southeast')

%% Frequency domain identification

nx_freq = 4;

OL_freqDat_trd      = idfrd(squeeze(sysID.nonPar.trd.OL.raw.ResponseData) ,sysID.nonPar.trd.OL.raw.Frequency ,1);
OL_freqDat_lpm      = idfrd(squeeze(sysID.nonPar.lpm.OL.raw.ResponseData) ,sysID.nonPar.lpm.OL.raw.Frequency ,1);
OL_freqDat_trd_filt = idfrd(squeeze(sysID.nonPar.trd.OL.filt.ResponseData),sysID.nonPar.trd.OL.filt.Frequency,1);
OL_freqDat_lpm_filt = idfrd(squeeze(sysID.nonPar.lpm.OL.filt.ResponseData),sysID.nonPar.lpm.OL.filt.Frequency,1);

% sysID.par.freqTRD.sim.OL.raw = ssest(OL_freqDat_trd,nx_freq,optSS_sim);
% sysID.par.freqLPM.sim.OL.raw = ssest(OL_freqDat_lpm,nx_freq,optSS_sim);
% sysID.par.freqTRD.sim.OL.raw  = ssest(OL_freqDat_trd     ,init_sys,optSS_sim);
% sysID.par.freqLPM.sim.OL.raw  = ssest(OL_freqDat_lpm     ,init_sys,optSS_sim);
sysID.par.freqTRD.sim.OL.raw  = ssest(OL_freqDat_trd     ,sysID.par.fixedOrder.sim.OL.raw,optSS_sim);
sysID.par.freqLPM.sim.OL.raw  = ssest(OL_freqDat_lpm     ,sysID.par.fixedOrder.sim.OL.raw,optSS_sim);
% sysID.par.freqTRD.sim.OL.filt = ssest(OL_freqDat_trd_filt,sysID.par.fixedOrder.sim.OL.raw,optSS_sim);
% sysID.par.freqLPM.sim.OL.filt = ssest(OL_freqDat_lpm_filt,sysID.par.fixedOrder.sim.OL.raw,optSS_sim);

%%
figure(baseFig_sim+701);clf
    % bode(sysID.par.freqTRD.sim.OL.raw,sysID.par.freqLPM.sim.OL.raw,sysID.par.freqTRD.sim.OL.filt,sysID.par.freqLPM.sim.OL.filt,G(1,1),'g--',opt);grid minor
    % bode(sysID.par.freqTRD.sim.OL.raw,sysID.par.freqLPM.sim.OL.raw,G(1,1),'g--',opt);grid minor
    % bode(sysID.par.freqTRD.sim.OL.raw,sysID.par.freqLPM.sim.OL.raw,sysID.nonPar.lpm.OL.raw,'g--',opt);grid minor
    bode(sysID.par.freqTRD.sim.OL.raw,sysID.par.freqLPM.sim.OL.raw,sysID.par.straight.sim.OL.raw,'g--',opt);grid minor

% figure(baseFig_sim+801);clf
%     resid(OL_freqDat_trd,sysID.par.freqTRD.sim.OL.raw,sysID.par.freqLPM.sim.OL.raw);grid minor
% figure(baseFig_sim+802);clf
%     resid(OL_freqDat_lpm,sysID.par.freqTRD.sim.OL.raw,sysID.par.freqLPM.sim.OL.raw);grid minor

% figure(baseFig_sim+801);clf
%     resid(OL_dat,sysID.par.freqTRD.sim.OL.raw,sysID.par.freqLPM.sim.OL.raw);grid minor
% figure(baseFig_sim+802);clf
%     resid(OL_dat,sysID.par.freqTRD.sim.OL.raw,sysID.par.freqLPM.sim.OL.raw);grid minor































