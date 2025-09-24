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
baseFig_sim = 8000;

makeValSet  = true;
removeTrans = false;

%% Make the system used for data generation
disp('Constructing system model')
main_lumpedSystem_water_V3

G = sysTMC(32,[3 6]);

%% Define some simulation variables
% fs = 1;     % [Hz]    Sampling frequency
fs = 0.01;     % [Hz]    Sampling frequency
Ts = 1/fs;  % [s]     Sampling time
tMeas = 8; % [hours] Measurement time

T0 = 23; % [degC] Inital Temperature
Tamb = 23; % [degC] Ambient Temperature

N = tMeas*3600*fs; % [-] Measurement samples
tVec = linspace(fs,tMeas*3600,N); % [s] Time vector

%% Define the noise on the signal
v     = 0.15.*randn(N,1);
v_val = 0.15.*randn(N,1);

%% Define the input signal
maxAmplitude = 24; % [W] Maximum exitation signal
positiveOnly = 1;  % [-] Define whether the exitation signal can be negative

disp('Generating multisine for identification')
dist     = genMultisine(fs, N, 1, maxAmplitude,     positiveOnly);
dist_val = dist;

%% Plot the noise signal
% Make the PSD of the noise and disturbance signals
    [amp_V,cas_V,freq_V] = fftCas_V2(v,fs);
    [amp_D,cas_D,freq_D] = fftCas_V2(dist,fs);

% Define FRD models of the PSD for easy use/plotting
    V      = frd(amp_V,freq_V);
    D      = frd(amp_D,freq_D);

% Visualize
    figure(baseFig_sim);clf
    set(gcf,'position',[700 100 700 500])
        subplot(211);
            semilogx(freq_V,db(amp_V),LineWidth=lw);grid minor
                xlabel('Frequency [Hz]')
                ylabel('Power [$degC^2/Hz$]')
                title('Power of noise signal added to simulation output')
        subplot(212)
            semilogx(freq_D,db(amp_D),LineWidth=lw);grid minor
                xlabel('Frequency [Hz]')
                ylabel('Power [$W^2/Hz$]')
                title('Power of disturbance/identification signal applied to simulation')

    figure(baseFig_sim+1);clf
    set(gcf,'position',[700 100 700 300])
        semilogx(freq_V,db(amp_V),LineWidth=lw);grid minor
            xlabel('Frequency [Hz]')
            ylabel('Power [$degC^2/Hz$]')
            title('Power of noise signal added to simulation output')

    figure(baseFig_sim+2);clf
    set(gcf,'position',[700 100 700 300])
        semilogx(freq_D,db(amp_D),LineWidth=lw);grid minor
            xlabel('Frequency [Hz]')
            ylabel('Power [$W^2/Hz$]')
            title('Power of disturbance/identification signal applied to simulation')

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

figure(baseFig_sim+101);clf
    set(gcf,'position',[700 100 700 300])
    plot(tVec,y_OL,LineWidth=lw);grid minor
        xlabel('Time [s]')
        ylabel('Temperature [degC]')
        title(['Thermal mass temperature, open loop identification, ',num2str(tMeas),' hours'])

%% Define the in and output vector sizes

% Define the range of the data that is used for identification and validation
if removeTrans
    idx1_OL = 8*3600*fs;
else
    idx1_OL = 1;
end
if makeValSet
    idx2_OL = N-(N-8*3600*fs)*0.2;
else
    idx2_OL = N;
end

idxRange_OL_trns = 1:idx1_OL+1;
idxRange_OL      = idx1_OL:idx2_OL;
% idxRange_OL_val  = (idx2_OL-1):N;
idxRange_OL_val  = (idx1_OL):N;

%% Define the in and output vectors, for easy and consistent use
% Open loop
outputData_OL_trns = y_OL(idxRange_OL_trns)-ambVec(idxRange_OL_trns)';
inputData_OL_trns  = dist(idxRange_OL_trns);

outputData_OL = y_OL(idxRange_OL)-ambVec(idxRange_OL)';
inputData_OL  = dist(idxRange_OL);

outputData_OL_val = y_OL_val(idxRange_OL_val)-ambVec(idxRange_OL_val)';
inputData_OL_val  = dist_val(idxRange_OL_val);

outputData_OL_full = y_OL-ambVec';
inputData_OL_full  = dist;

% Define the time vectors
tVec_OL_trns = tVec(idxRange_OL_trns);
tVec_OL      = tVec(idxRange_OL);
tVec_OL_val  = tVec(idxRange_OL_val);
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
    
%% FRF measurements using the new state vector
% Make the FRF's using the raw data
    [sysID.nonPar.trd.OL.raw, ~] = makeOpenLoopFRF_sysIdent(outputData_OL, inputData_OL(1,:), fs);

% Make the FRF's using the filtered data
    [sysID.nonPar.trd.OL.filt, ~] = makeOpenLoopFRF_sysIdent(outputData_OL_filt, inputData_OL_filt(1,:), fs);

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
    set(gcf,'position',[500 100 900 500])
    sgtitle(['Heater to TM temperature, open loop identification using unfiltered and filtered data, ',num2str(tMeas),' hours'])
        subplot(221)
            simpleBodemag(sysID.nonPar.trd.OL.raw,'Hz',lw);hold on;grid minor
            simpleBodemag(sysID.nonPar.lpm.OL.raw,'Hz',lw);grid minor
            simpleBodemag(G(1,1),'Hz',lw,'g--');grid minor
            xlim([1e-4 0.01])
                title(['Using unfiltered data'])
                legend('Traditional','LPM','Model','location','best')
        subplot(223)
            simpleBodephase(sysID.nonPar.trd.OL.raw,'Hz',lw,'wrap');hold on;grid minor
            simpleBodephase(sysID.nonPar.lpm.OL.raw,'Hz',lw,'wrap');grid minor
            simpleBodephase(G(1,1),'Hz',lw,'g--','wrap');grid minor
            xlim([1e-4 0.01])
                legend('Traditional','LPM','Model','location','best')
        subplot(222)
            simpleBodemag(sysID.nonPar.trd.OL.filt,'Hz',lw);hold on;grid minor
            simpleBodemag(sysID.nonPar.lpm.OL.filt,'Hz',lw);grid minor
            simpleBodemag(G(1,1),'Hz',lw,'g--');grid minor
            xlim([1e-4 0.01])
                title(['Using filtered data'])
                legend('Traditional','LPM','Model','location','best')
        subplot(224)
            simpleBodephase(sysID.nonPar.trd.OL.filt,'Hz',lw,'wrap');hold on;grid minor
            simpleBodephase(sysID.nonPar.lpm.OL.filt,'Hz',lw,'wrap');grid minor
            simpleBodephase(G(1,1),'Hz',lw,'g--','wrap');grid minor
            xlim([1e-4 0.01])
                legend('Traditional','LPM','Model','location','best')

%% First order parametric approximation, using time data
    sysID.par.firstAprrox.OL.raw  = step_sysID(inputData_OL_full',zeros(size(inputData_OL_full))',outputData_OL_full,tVec_OL_full,10000*fs);
    sysID.par.firstAprrox.OL.filt = step_sysID(inputData_OL_full_filt',zeros(size(inputData_OL_full_filt))',outputData_OL_full_filt',tVec_OL_full,10000*fs);

%% Visualization of the first order approximation
bodeRange = logspace(-6,-1,100);

figure(baseFig_sim+203);clf
    set(gcf,'position',[500 100 900 500])
    sgtitle(['Heater to thermal mass temperature, first order approximation, ',num2str(tMeas),' hours'])
        subplot(211)
            simpleBodemag(sysID.par.firstAprrox.OL.raw,'Hz',lw ,bodeRange);hold on;grid minor
            simpleBodemag(sysID.par.firstAprrox.OL.filt,'Hz',lw,bodeRange,'-.');grid minor
            simpleBodemag(G(1,1),'Hz',lw,'g--',bodeRange);grid minor
                legend('Gotten from unfiltered data','Gotten from filtered data','Model','location','best')
        subplot(212)
            simpleBodephase(sysID.par.firstAprrox.OL.raw,'Hz',lw,'wrap' ,bodeRange);hold on;grid minor
            simpleBodephase(sysID.par.firstAprrox.OL.filt,'Hz',lw,'wrap',bodeRange,'-.');grid minor
            simpleBodephase(G(1,1),'Hz',lw,'g--','wrap',bodeRange);grid minor
                legend('Gotten from unfiltered data','Gotten from filtered data','Model','location','best')

%% Parametric system identification, using time data, pre-requisites

% Definitions for identification
    nx = 3; % Model order for fixed order identification

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

% All tested parametric identification options using OL data, using prediction focus
        disp('Open Loop identification using an initial system: in progress')
    sysID.par.initSys.sim.OL.raw  = ssest(OL_dat,init_sys,optSS_sim);
    sysID.par.initSys.sim.OL.filt = ssest(OL_dat_filt,init_sys,optSS_sim);
        disp('Open Loop identification using an initial system: done')

        disp(['Open Loop identification using a fixed order (nx=',num2str(nx),'): in progress'])
    sysID.par.fixedOrder.sim.OL.raw  = ssest(OL_dat,nx,optSS_sim);
    sysID.par.fixedOrder.sim.OL.filt = ssest(OL_dat_filt,nx,optSS_sim);
        disp(['Open Loop identification using a fixed order (nx=',num2str(nx),'): done'])

% All tested parametric identification options using OL data, using prediction focus
        disp('Open Loop identification using an initial system: in progress')
    sysID.par.initSys.prd.OL.raw  = ssest(OL_dat,init_sys,optSS_prd);
    sysID.par.initSys.prd.OL.filt = ssest(OL_dat_filt,init_sys,optSS_prd);
        disp('Open Loop identification using an initial system: done')

        disp(['Open Loop identification using a fixed order (nx=',num2str(nx),'): in progress'])
    sysID.par.fixedOrder.prd.OL.raw  = ssest(OL_dat,nx,optSS_prd);
    sysID.par.fixedOrder.prd.OL.filt = ssest(OL_dat_filt,nx,optSS_prd);
        disp(['Open Loop identification using a fixed order (nx=',num2str(nx),'): done'])


        disp('Open Loop identification using the raw data: in progress')
    sysID.par.straight.sim.OL.raw  = ssest(OL_dat);
    sysID.par.straight.sim.OL.filt = ssest(OL_dat_filt);
        disp('Open Loop identification using the raw data: done')

%% Visualize and compare the identifications in frequency domain
bodeRange = logspace(-6,-1,1000);

figure(baseFig_sim+301);clf
    set(gcf,'position',[500 100 700 450])
    sgtitle('Bode plots of identification using raw Open Loop data')
        subplot(211)
            simpleBodemag(sysID.par.initSys.sim.OL.raw   ,'Hz',lw,bodeRange);hold on;grid minor
            simpleBodemag(sysID.par.fixedOrder.sim.OL.raw,'Hz',lw,bodeRange);grid minor
            simpleBodemag(sysID.par.firstAprrox.OL.raw   ,'Hz',lw,bodeRange);grid minor
            simpleBodemag(sysID.par.straight.sim.OL.raw  ,'Hz',lw,bodeRange);grid minor
            simpleBodemag(sysTMC(32,3)                   ,'Hz',lw,'g--',bodeRange);grid minor
                legend('InitSysIdent','Fixed order','First order approx','Straight data','Model','location','best')
        subplot(212)
            simpleBodephase(sysID.par.initSys.sim.OL.raw   ,'Hz',lw,bodeRange,'wrap');hold on;grid minor
            simpleBodephase(sysID.par.fixedOrder.sim.OL.raw,'Hz',lw,bodeRange,'wrap');grid minor
            simpleBodephase(sysID.par.firstAprrox.OL.raw   ,'Hz',lw,bodeRange,'wrap');grid minor
            simpleBodephase(sysID.par.straight.sim.OL.raw  ,'Hz',lw,bodeRange,'wrap');grid minor
            simpleBodephase(sysTMC(32,3)                   ,'Hz',lw,'g--',bodeRange,'wrap');grid minor
                legend('InitSysIdent','Fixed order','First order approx','Straight data','Model','location','best')

figure(baseFig_sim+302);clf
    set(gcf,'position',[500 100 700 450])
    sgtitle('Bode plots of identification using filtered Open Loop data')
        subplot(211)
            simpleBodemag(sysID.par.initSys.sim.OL.filt   ,'Hz',lw,bodeRange);hold on;grid minor
            simpleBodemag(sysID.par.fixedOrder.sim.OL.filt,'Hz',lw,bodeRange);grid minor
            simpleBodemag(sysID.par.firstAprrox.OL.filt   ,'Hz',lw,bodeRange);grid minor
            simpleBodemag(sysID.par.straight.sim.OL.filt  ,'Hz',lw,bodeRange);grid minor
            simpleBodemag(sysTMC(32,3)                    ,'Hz',lw,'g--',bodeRange);grid minor
                legend('InitSysIdent','Fixed order','First order approx','Straight data','Model','location','best')
        subplot(212)
            simpleBodephase(sysID.par.initSys.sim.OL.filt   ,'Hz',lw,bodeRange,'wrap');hold on;grid minor
            simpleBodephase(sysID.par.fixedOrder.sim.OL.filt,'Hz',lw,bodeRange,'wrap');grid minor
            simpleBodephase(sysID.par.firstAprrox.OL.filt   ,'Hz',lw,bodeRange,'wrap');grid minor
            simpleBodephase(sysID.par.straight.sim.OL.filt  ,'Hz',lw,bodeRange,'wrap');grid minor
            simpleBodephase(sysTMC(32,3)                    ,'Hz',lw,'g--',bodeRange,'wrap');grid minor
                legend('InitSysIdent','Fixed order','First order approx','Straight data','Model','location','best')

% figure(baseFig_sim+303);clf
%     set(gcf,'position',[100 60 1200 800])
%     sgtitle('Bode plots of identification using filtered Open Loop data')
%         subplot(421)
%             simpleBodemag(sysID.par.initSys.sim.OL.raw ,'Hz',lw,bodeRange);hold on;grid minor
%             simpleBodemag(sysID.par.initSys.sim.OL.filt,'Hz',lw,bodeRange);grid minor
%             simpleBodemag(sysTMC(32,3)                 ,'Hz',lw,'g--',bodeRange);grid minor
%             xlim([bodeRange(1) bodeRange(end)])
%                 title('Initial system')
%                 legend('Raw','Filtered','location','best')
%         subplot(423)
%             simpleBodephase(sysID.par.initSys.sim.OL.raw ,'Hz',lw,bodeRange,'wrap');hold on;grid minor
%             simpleBodephase(sysID.par.initSys.sim.OL.filt,'Hz',lw,bodeRange,'wrap');grid minor
%             simpleBodephase(sysTMC(32,3)                 ,'Hz',lw,'g--',bodeRange,'wrap');grid minor
%             xlim([bodeRange(1) bodeRange(end)])
%                 legend('Raw','Filtered','location','best')
%         subplot(422)
%             simpleBodemag(sysID.par.fixedOrder.sim.OL.raw ,'Hz',lw,bodeRange);hold on;grid minor
%             simpleBodemag(sysID.par.fixedOrder.sim.OL.filt,'Hz',lw,bodeRange);grid minor
%             simpleBodemag(sysTMC(32,3)                    ,'Hz',lw,'g--',bodeRange);grid minor
%             xlim([bodeRange(1) bodeRange(end)])
%                 title(['Fixed order, nx=',num2str(nx)])
%                 legend('Raw','Filtered','location','best')
%         subplot(424)
%             simpleBodephase(sysID.par.fixedOrder.sim.OL.raw ,'Hz',lw,bodeRange,'wrap');hold on;grid minor
%             simpleBodephase(sysID.par.fixedOrder.sim.OL.filt,'Hz',lw,bodeRange,'wrap');grid minor
%             simpleBodephase(sysTMC(32,3)                    ,'Hz',lw,'g--',bodeRange,'wrap');grid minor
%             xlim([bodeRange(1) bodeRange(end)])
%                 legend('Raw','Filtered','location','best')
%         subplot(425)
%             simpleBodemag(sysID.par.firstAprrox.OL.raw ,'Hz',lw,bodeRange);hold on;grid minor
%             simpleBodemag(sysID.par.firstAprrox.OL.filt,'Hz',lw,bodeRange);grid minor
%             simpleBodemag(sysTMC(32,3)                 ,'Hz',lw,'g--',bodeRange);grid minor
%             xlim([bodeRange(1) bodeRange(end)])
%                 title('First order approximation')
%                 legend('Raw','Filtered','location','best')
%         subplot(427)
%             simpleBodephase(sysID.par.firstAprrox.OL.raw ,'Hz',lw,bodeRange,'wrap');hold on;grid minor
%             simpleBodephase(sysID.par.firstAprrox.OL.filt,'Hz',lw,bodeRange,'wrap');grid minor
%             simpleBodephase(sysTMC(32,3)                 ,'Hz',lw,'g--',bodeRange,'wrap');grid minor
%             xlim([bodeRange(1) bodeRange(end)])
%                 legend('Raw','Filtered','location','best')
%         subplot(426)
%             simpleBodemag(sysID.par.straight.sim.OL.raw ,'Hz',lw,bodeRange);hold on;grid minor
%             simpleBodemag(sysID.par.straight.sim.OL.filt,'Hz',lw,bodeRange);grid minor
%             simpleBodemag(sysTMC(32,3)                  ,'Hz',lw,'g--',bodeRange);grid minor
%             xlim([bodeRange(1) bodeRange(end)])
%                 title('Straight data')
%                 legend('Raw','Filtered','location','best')
%         subplot(428)
%             simpleBodephase(sysID.par.straight.sim.OL.raw ,'Hz',lw,bodeRange,'wrap');hold on;grid minor
%             simpleBodephase(sysID.par.straight.sim.OL.filt,'Hz',lw,bodeRange,'wrap');grid minor
%             simpleBodephase(sysTMC(32,3)                  ,'Hz',lw,'g--',bodeRange,'wrap');grid minor
%             xlim([bodeRange(1) bodeRange(end)])
%                 legend('Raw','Filtered','location','best')

figure(baseFig_sim+304);clf
    set(gcf,'position',[100 60 1400 450])
    sgtitle('Bode plots of identification using filtered Open Loop data')
        subplot(241)
            simpleBodemag(sysID.par.initSys.sim.OL.raw ,'Hz',lw,bodeRange);hold on;grid minor
            simpleBodemag(sysID.par.initSys.sim.OL.filt,'Hz',lw,bodeRange);grid minor
            simpleBodemag(sysTMC(32,3)                 ,'Hz',lw,'g--',bodeRange);grid minor
            xlim([bodeRange(1) bodeRange(end)])
                title('Initial system')
                legend('Using raw data','Using filtered data','location','best')
        subplot(245)
            simpleBodephase(sysID.par.initSys.sim.OL.raw ,'Hz',lw,bodeRange,'wrap');hold on;grid minor
            simpleBodephase(sysID.par.initSys.sim.OL.filt,'Hz',lw,bodeRange,'wrap');grid minor
            simpleBodephase(sysTMC(32,3)                 ,'Hz',lw,'g--',bodeRange,'wrap');grid minor
            xlim([bodeRange(1) bodeRange(end)])
                legend('Using raw data','Using filtered data','location','best')
        subplot(242)
            simpleBodemag(sysID.par.fixedOrder.sim.OL.raw ,'Hz',lw,bodeRange);hold on;grid minor
            simpleBodemag(sysID.par.fixedOrder.sim.OL.filt,'Hz',lw,bodeRange);grid minor
            simpleBodemag(sysTMC(32,3)                    ,'Hz',lw,'g--',bodeRange);grid minor
            xlim([bodeRange(1) bodeRange(end)])
                title(['Fixed order, nx=',num2str(nx)])
                legend('Using raw data','Using filtered data','location','best')
        subplot(246)
            simpleBodephase(sysID.par.fixedOrder.sim.OL.raw ,'Hz',lw,bodeRange,'wrap');hold on;grid minor
            simpleBodephase(sysID.par.fixedOrder.sim.OL.filt,'Hz',lw,bodeRange,'wrap');grid minor
            simpleBodephase(sysTMC(32,3)                    ,'Hz',lw,'g--',bodeRange,'wrap');grid minor
            xlim([bodeRange(1) bodeRange(end)])
                legend('Using raw data','Using filtered data','location','best')
        subplot(243)
            simpleBodemag(sysID.par.firstAprrox.OL.raw ,'Hz',lw,bodeRange);hold on;grid minor
            simpleBodemag(sysID.par.firstAprrox.OL.filt,'Hz',lw,bodeRange);grid minor
            simpleBodemag(sysTMC(32,3)                 ,'Hz',lw,'g--',bodeRange);grid minor
            xlim([bodeRange(1) bodeRange(end)])
                title('First order approximation')
                legend('Using raw data','Using filtered data','location','best')
        subplot(247)
            simpleBodephase(sysID.par.firstAprrox.OL.raw ,'Hz',lw,bodeRange,'wrap');hold on;grid minor
            simpleBodephase(sysID.par.firstAprrox.OL.filt,'Hz',lw,bodeRange,'wrap');grid minor
            simpleBodephase(sysTMC(32,3)                 ,'Hz',lw,'g--',bodeRange,'wrap');grid minor
            xlim([bodeRange(1) bodeRange(end)])
                legend('Using raw data','Using filtered data','location','best')
        subplot(244)
            simpleBodemag(sysID.par.straight.sim.OL.raw ,'Hz',lw,bodeRange);hold on;grid minor
            simpleBodemag(sysID.par.straight.sim.OL.filt,'Hz',lw,bodeRange);grid minor
            simpleBodemag(sysTMC(32,3)                  ,'Hz',lw,'g--',bodeRange);grid minor
            xlim([bodeRange(1) bodeRange(end)])
                title('Straight data')
                legend('Using raw data','Using filtered data','location','best')
        subplot(248)
            simpleBodephase(sysID.par.straight.sim.OL.raw ,'Hz',lw,bodeRange,'wrap');hold on;grid minor
            simpleBodephase(sysID.par.straight.sim.OL.filt,'Hz',lw,bodeRange,'wrap');grid minor
            simpleBodephase(sysTMC(32,3)                  ,'Hz',lw,'g--',bodeRange,'wrap');grid minor
            xlim([bodeRange(1) bodeRange(end)])
                legend('Using raw data','Using filtered data','location','best')

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
    OLsim_sysID_predTests_V3

% Run the residual tests and plot them
    OLsim_sysID_resTests_V3


%% Frequency domain identification

nx_freq = 4;

OL_freqDat_trd      = idfrd(squeeze(sysID.nonPar.trd.OL.raw.ResponseData) ,sysID.nonPar.trd.OL.raw.Frequency ,1);
OL_freqDat_lpm      = idfrd(squeeze(sysID.nonPar.lpm.OL.raw.ResponseData) ,sysID.nonPar.lpm.OL.raw.Frequency ,1);
OL_freqDat_trd_filt = idfrd(squeeze(sysID.nonPar.trd.OL.filt.ResponseData),sysID.nonPar.trd.OL.filt.Frequency,1);
OL_freqDat_lpm_filt = idfrd(squeeze(sysID.nonPar.lpm.OL.filt.ResponseData),sysID.nonPar.lpm.OL.filt.Frequency,1);

% sysID.par.freqTRD.sim.OL.raw = ssest(OL_freqDat_trd,nx_freq,optSS_sim);
% sysID.par.freqLPM.sim.OL.raw = ssest(OL_freqDat_lpm,nx_freq,optSS_sim);
sysID.par.freqTRD.sim.OL.raw  = ssest(OL_freqDat_trd     ,init_sys,optSS_sim);
sysID.par.freqLPM.sim.OL.raw  = ssest(OL_freqDat_lpm     ,init_sys,optSS_sim);
% sysID.par.freqTRD.sim.OL.raw  = ssest(OL_freqDat_trd     ,sysID.par.fixedOrder.sim.OL.raw,optSS_sim);
% sysID.par.freqLPM.sim.OL.raw  = ssest(OL_freqDat_lpm     ,sysID.par.fixedOrder.sim.OL.raw,optSS_sim);
% sysID.par.freqTRD.sim.OL.filt = ssest(OL_freqDat_trd_filt,sysID.par.fixedOrder.sim.OL.raw,optSS_sim);
% sysID.par.freqLPM.sim.OL.filt = ssest(OL_freqDat_lpm_filt,sysID.par.fixedOrder.sim.OL.raw,optSS_sim);

%%
figure(baseFig_sim+701);clf
    % bode(sysID.par.freqTRD.sim.OL.raw,sysID.par.freqLPM.sim.OL.raw,sysID.par.freqTRD.sim.OL.filt,sysID.par.freqLPM.sim.OL.filt,G(1,1),'g--',opt);grid minor
    bode(sysID.par.freqTRD.sim.OL.raw,sysID.par.freqLPM.sim.OL.raw,G(1,1),'g--',opt);grid minor






























