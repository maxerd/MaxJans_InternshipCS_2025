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
Tamb = 23; % [degC] Ambient Temperature

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


%% Some user inputs
makeValSet  = true;  % [true/false] Is a seperate validation set needed?
removeTrans = false; % [true/false] Does the transient need to be removed

%% Define some simulation variables

T0 = 23; % [degC] Inital Temperature

N = tMeas*3600*fs; % [-] Measurement samples
tVec = linspace(fs,tMeas*3600,N); % [s] Time vector

%% Define the noise on the signal
v     = 0.15.*randn(N,1);
v_val = 0.15.*randn(N,1);

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

%% First order parametric approximation, using time data
    sysID.par.firstAprrox.OL.raw = step_sysID(dist',zeros(size(dist))',y_OL-ambVec'',tVec,10000*fs);
    sysID.par.firstAprrox.OL.filt = step_sysID(dist',zeros(size(dist))',y_OL-ambVec'',tVec,10000*fs);

%% Simulate the systems

disp('Running simulation')

% Ambient temperature vector, constant ambient temperature for simulations
    ambVec = Tamb.*ones(1,length(dist));

% Define the initial temperatures (states) for the simulation, also constant.
    x0_OL = T0.*ones(size(G.A,1),1);

% Run the open loop simulation, with the previously defined noise signal added
    y_OL     = lsim(G,[dist    ;ambVec],tVec,x0_OL)+v;
    y_OL_val = lsim(G,[dist_val;ambVec],tVec,x0_OL)+v_val;

% First order parametric approximation, using time data
    sysID.par.firstAprrox.OL.raw = step_sysID(dist',zeros(size(dist))',y_OL-ambVec'',tVec,10000*fs);

% Remove the transient using the first order approximation
    y_approx1     = lsim(sysID.par.firstAprrox.OL.raw,[dist],tVec);
    y_approx1_val = lsim(sysID.par.firstAprrox.OL.raw,[dist_val],tVec);
    % figure;plot(tVec,y_OL-y_approx1)
    % figure;plot(tVec,y_OL,tVec,y_approx1)
    % y_OL     = ((y_OL-y_approx1));
    figure;plot(tVec,y_OL);grid minor
    % y_OL_val = (y_OL_val-y_approx1_val);

%% Plot the simulation results

figure(baseFig_sim+101);clf
    % set(gcf,'position',[700 100 700 300])
    plot(tVec,y_OL-detrend(y_OL,6),LineWidth=lw);grid minor;hold on
    % plot(tVec,v);grid minor
    % xlim([0 1000])
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
    idx2_OL = N-(N-6*3600*fs)*0.2;
else
    idx2_OL = N;
end

idxRange_OL_trns = 1:idx1_OL+1;
idxRange_OL      = idx1_OL:idx2_OL;
% idxRange_OL_val  = (idx2_OL-1):N;
idxRange_OL_val  = (idx2_OL):N;

%% Define the in and output vectors, for easy and consistent use
% Open loop
sysID.data.OL.trans.out = y_OL(idxRange_OL_trns)-ambVec(idxRange_OL_trns)';
sysID.data.OL.trans.in  = dist(idxRange_OL_trns);

sysID.data.OL.train.out = y_OL(idxRange_OL)-ambVec(idxRange_OL)';
sysID.data.OL.train.in  = dist(idxRange_OL);

sysID.data.OL.val.out   = y_OL_val(idxRange_OL_val)-ambVec(idxRange_OL_val)';
sysID.data.OL.val.in    = dist_val(idxRange_OL_val);

sysID.data.OL.full.out  = y_OL-ambVec';
sysID.data.OL.full.in   = dist;

% Define the time vectors
sysID.data.OL.trans.tVec = tVec(idxRange_OL_trns);
sysID.data.OL.train.tVec = tVec(idxRange_OL);
sysID.data.OL.val.tVec   = tVec(idxRange_OL_val);
sysID.data.OL.full.tVec  = tVec;

% tVec_CL      = tVec(idxRange_CL);

%% Filter the data vectors to remove some of the effect of noise
% High order lowpass filter, to filter away all the noise above the 
% multisine exitation frequency
    LPfilt = orderLP(0.01,2,0,0,0.7)^2;
    % [LPfiltB,LPfiltA] = butter(4,0.005);
    % LPfilt = tf(LPfiltB,LPfiltA,1);

% Actually filter the data
    % Open loop data
    sysID.dataFilt.OL.trans.out = lsim(LPfilt,sysID.data.OL.trans.out-sysID.data.OL.trans.out(1),sysID.data.OL.trans.tVec)'+sysID.data.OL.trans.out(1);
    sysID.dataFilt.OL.trans.in = lsim(LPfilt,sysID.data.OL.trans.in-sysID.data.OL.trans.in(1),sysID.data.OL.trans.tVec)'+sysID.data.OL.trans.in(1);

    sysID.dataFilt.OL.train.out = lsim(LPfilt,sysID.data.OL.train.out-sysID.data.OL.train.out(1),sysID.data.OL.train.tVec)'+sysID.data.OL.train.out(1);
    sysID.dataFilt.OL.train.in = lsim(LPfilt,sysID.data.OL.train.in-sysID.data.OL.train.in(1),sysID.data.OL.train.tVec)'+sysID.data.OL.train.in(1);

    sysID.dataFilt.OL.val.out = lsim(LPfilt,sysID.data.OL.val.out-sysID.data.OL.val.out(1),sysID.data.OL.val.tVec)'+sysID.data.OL.val.out(1);
    sysID.dataFilt.OL.val.in = lsim(LPfilt,sysID.data.OL.val.in-sysID.data.OL.val.in(1),sysID.data.OL.val.tVec)'+sysID.data.OL.val.in(1);

    sysID.dataFilt.OL.full.out = lsim(LPfilt,sysID.data.OL.full.out-sysID.data.OL.full.out(1),sysID.data.OL.full.tVec)'+sysID.data.OL.full.out(1);
    sysID.dataFilt.OL.full.in = lsim(LPfilt,sysID.data.OL.full.in-sysID.data.OL.full.in(1),sysID.data.OL.full.tVec)'+sysID.data.OL.full.in(1);
    
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
    print('./_Figures/OL/nonPar/tradAndLPM_FRF_1','-depsc')

%% First order parametric approximation, using time data
    % sysID.par.firstAprrox.OL.raw  = step_sysID(sysID.data.OL.full.in',zeros(size(sysID.data.OL.full.in))',sysID.data.OL.full.out,sysID.data.OL.full.tVec,10000*fs);
    % sysID.par.firstAprrox.OL.filt = step_sysID(sysID.dataFilt.OL.full.in',zeros(size(sysID.dataFilt.OL.full.in))',sysID.dataFilt.OL.full.out',sysID.data.OL.full.tVec,10000*fs);
    sysID.par.firstAprrox.OL.raw  = firstOrderApprox(sysID.data.OL.train.id);
    sysID.par.firstAprrox.OL.filt = firstOrderApprox(sysID.dataFilt.OL.train.id);

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
    print('./_Figures/OL/parTime/firstOrderApprox_1','-depsc')

%% Parametric system identification, using time data, pre-requisites

% Definitions for identification
    nx = 5; % Model order for fixed order identification

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
    sysID.data.OL.train.id      = iddata(sysID.data.OL.train.out          ,sysID.data.OL.train.in(:,:)'    ,Ts);
    sysID.dataFilt.OL.train.id = iddata(sysID.dataFilt.OL.train.out'    ,sysID.dataFilt.OL.train.in'    ,Ts);

% Define the data that is not used in the identification as validation data

if makeValSet
    sysID.data.OL.val.id      = iddata(sysID.data.OL.val.out      ,sysID.data.OL.val.in',Ts);
    sysID.dataFilt.OL.val.id = iddata(sysID.dataFilt.OL.val.out',sysID.data.OL.val.in',Ts);
else
    sysID.data.OL.val.id      = iddata(1,1,Ts);
    sysID.dataFilt.OL.val.id = iddata(1,1,Ts);
end
if removeTrans
    sysID.data.OL.trans.id      = iddata(sysID.data.OL.trans.out      ,sysID.data.OL.trans.in',Ts);
    sysID.dataFilt.OL.trans.id = iddata(sysID.dataFilt.OL.trans.out',sysID.data.OL.trans.in',Ts);
else
    sysID.data.OL.trans.id      = iddata(1,1,Ts);
    sysID.dataFilt.OL.trans.id = iddata(1,1,Ts);
end

optSS_sim = ssestOptions('Focus','Simulation');
optSS_prd = ssestOptions('Focus','Prediction');

%% Parametric system identification, using time data

% opts.focus = 'sim';
% opts.bodeRange = logspace(-6,-1,1000);
% opts.lw = 1.5;
figure;initializedIdent_OL(sysID.data.OL.train.id,sysID.par.firstAprrox.OL.raw,opts);hold on
        fixedOrderIdent_OL(sysID.data.OL.train.id,3,opts)
               rawIdent_OL(sysID.data.OL.train.id,opts)
%%

% All tested parametric identification options using OL data, using prediction focus
        disp('Open Loop identification using an initial system: in progress')
    sysID.par.initSys.sim.OL.raw  = ssest(sysID.data.OL.train.id,init_sys,optSS_sim);
    sysID.par.initSys.sim.OL.filt = ssest(sysID.dataFilt.OL.train.id,init_sys,optSS_sim);
        disp('Open Loop identification using an initial system: done')

        disp(['Open Loop identification using a fixed order (nx=',num2str(nx),'): in progress'])
    sysID.par.fixedOrder.sim.OL.raw  = ssest(sysID.data.OL.train.id,nx,optSS_sim);
    sysID.par.fixedOrder.sim.OL.filt = ssest(sysID.dataFilt.OL.train.id,nx,optSS_sim);
        disp(['Open Loop identification using a fixed order (nx=',num2str(nx),'): done'])

% All tested parametric identification options using OL data, using prediction focus
        disp('Open Loop identification using an initial system: in progress')
    sysID.par.initSys.prd.OL.raw  = ssest(sysID.data.OL.train.id,init_sys,optSS_prd);
    sysID.par.initSys.prd.OL.filt = ssest(sysID.dataFilt.OL.train.id,init_sys,optSS_prd);
        disp('Open Loop identification using an initial system: done')

        disp(['Open Loop identification using a fixed order (nx=',num2str(nx),'): in progress'])
    sysID.par.fixedOrder.prd.OL.raw  = ssest(sysID.data.OL.train.id,nx,optSS_prd);
    sysID.par.fixedOrder.prd.OL.filt = ssest(sysID.dataFilt.OL.train.id,nx,optSS_prd);
        disp(['Open Loop identification using a fixed order (nx=',num2str(nx),'): done'])


        disp('Open Loop identification using the raw data: in progress')
    sysID.par.straight.sim.OL.raw  = ssest(sysID.data.OL.train.id);
    sysID.par.straight.sim.OL.filt = ssest(sysID.dataFilt.OL.train.id);
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
            simpleBodemag(G(1,1)                   ,'Hz',lw,'g--',bodeRange);grid minor
                legend('InitSysIdent','Fixed order','First order approx','Straight data','Model','location','best')
        subplot(212)
            simpleBodephase(sysID.par.initSys.sim.OL.raw   ,'Hz',lw,bodeRange,'wrap');hold on;grid minor
            simpleBodephase(sysID.par.fixedOrder.sim.OL.raw,'Hz',lw,bodeRange,'wrap');grid minor
            simpleBodephase(sysID.par.firstAprrox.OL.raw   ,'Hz',lw,bodeRange,'wrap');grid minor
            simpleBodephase(sysID.par.straight.sim.OL.raw  ,'Hz',lw,bodeRange,'wrap');grid minor
            simpleBodephase(G(1,1)                   ,'Hz',lw,'g--',bodeRange,'wrap');grid minor
                legend('InitSysIdent','Fixed order','First order approx','Straight data','Model','location','best')
    print('./_Figures/OL/parTime/compareBodeAll_raw_1','-depsc')

figure(baseFig_sim+302);clf
    set(gcf,'position',[500 100 700 450])
    sgtitle('Bode plots of identification using filtered Open Loop data')
        subplot(211)
            simpleBodemag(sysID.par.initSys.sim.OL.filt   ,'Hz',lw,bodeRange);hold on;grid minor
            simpleBodemag(sysID.par.fixedOrder.sim.OL.filt,'Hz',lw,bodeRange);grid minor
            simpleBodemag(sysID.par.firstAprrox.OL.filt   ,'Hz',lw,bodeRange);grid minor
            simpleBodemag(sysID.par.straight.sim.OL.filt  ,'Hz',lw,bodeRange);grid minor
            simpleBodemag(G(1,1)                    ,'Hz',lw,'g--',bodeRange);grid minor
                legend('InitSysIdent','Fixed order','First order approx','Straight data','Model','location','best')
        subplot(212)
            simpleBodephase(sysID.par.initSys.sim.OL.filt   ,'Hz',lw,bodeRange,'wrap');hold on;grid minor
            simpleBodephase(sysID.par.fixedOrder.sim.OL.filt,'Hz',lw,bodeRange,'wrap');grid minor
            simpleBodephase(sysID.par.firstAprrox.OL.filt   ,'Hz',lw,bodeRange,'wrap');grid minor
            simpleBodephase(sysID.par.straight.sim.OL.filt  ,'Hz',lw,bodeRange,'wrap');grid minor
            simpleBodephase(G(1,1)                    ,'Hz',lw,'g--',bodeRange,'wrap');grid minor
                legend('InitSysIdent','Fixed order','First order approx','Straight data','Model','location','best')
    print('./_Figures/OL/parTime/compareBodeAll_filt_1','-depsc')

% figure(baseFig_sim+303);clf
%     set(gcf,'position',[100 60 1200 800])
%     sgtitle('Bode plots of identification using filtered Open Loop data')
%         subplot(421)
%             simpleBodemag(sysID.par.initSys.sim.OL.raw ,'Hz',lw,bodeRange);hold on;grid minor
%             simpleBodemag(sysID.par.initSys.sim.OL.filt,'Hz',lw,bodeRange);grid minor
%             simpleBodemag(G(1,1)                 ,'Hz',lw,'g--',bodeRange);grid minor
%             xlim([bodeRange(1) bodeRange(end)])
%                 title('Initial system')
%                 legend('Raw','Filtered','location','best')
%         subplot(423)
%             simpleBodephase(sysID.par.initSys.sim.OL.raw ,'Hz',lw,bodeRange,'wrap');hold on;grid minor
%             simpleBodephase(sysID.par.initSys.sim.OL.filt,'Hz',lw,bodeRange,'wrap');grid minor
%             simpleBodephase(G(1,1)                 ,'Hz',lw,'g--',bodeRange,'wrap');grid minor
%             xlim([bodeRange(1) bodeRange(end)])
%                 legend('Raw','Filtered','location','best')
%         subplot(422)
%             simpleBodemag(sysID.par.fixedOrder.sim.OL.raw ,'Hz',lw,bodeRange);hold on;grid minor
%             simpleBodemag(sysID.par.fixedOrder.sim.OL.filt,'Hz',lw,bodeRange);grid minor
%             simpleBodemag(G(1,1)                    ,'Hz',lw,'g--',bodeRange);grid minor
%             xlim([bodeRange(1) bodeRange(end)])
%                 title(['Fixed order, nx=',num2str(nx)])
%                 legend('Raw','Filtered','location','best')
%         subplot(424)
%             simpleBodephase(sysID.par.fixedOrder.sim.OL.raw ,'Hz',lw,bodeRange,'wrap');hold on;grid minor
%             simpleBodephase(sysID.par.fixedOrder.sim.OL.filt,'Hz',lw,bodeRange,'wrap');grid minor
%             simpleBodephase(G(1,1)                    ,'Hz',lw,'g--',bodeRange,'wrap');grid minor
%             xlim([bodeRange(1) bodeRange(end)])
%                 legend('Raw','Filtered','location','best')
%         subplot(425)
%             simpleBodemag(sysID.par.firstAprrox.OL.raw ,'Hz',lw,bodeRange);hold on;grid minor
%             simpleBodemag(sysID.par.firstAprrox.OL.filt,'Hz',lw,bodeRange);grid minor
%             simpleBodemag(G(1,1)                 ,'Hz',lw,'g--',bodeRange);grid minor
%             xlim([bodeRange(1) bodeRange(end)])
%                 title('First order approximation')
%                 legend('Raw','Filtered','location','best')
%         subplot(427)
%             simpleBodephase(sysID.par.firstAprrox.OL.raw ,'Hz',lw,bodeRange,'wrap');hold on;grid minor
%             simpleBodephase(sysID.par.firstAprrox.OL.filt,'Hz',lw,bodeRange,'wrap');grid minor
%             simpleBodephase(G(1,1)                 ,'Hz',lw,'g--',bodeRange,'wrap');grid minor
%             xlim([bodeRange(1) bodeRange(end)])
%                 legend('Raw','Filtered','location','best')
%         subplot(426)
%             simpleBodemag(sysID.par.straight.sim.OL.raw ,'Hz',lw,bodeRange);hold on;grid minor
%             simpleBodemag(sysID.par.straight.sim.OL.filt,'Hz',lw,bodeRange);grid minor
%             simpleBodemag(G(1,1)                  ,'Hz',lw,'g--',bodeRange);grid minor
%             xlim([bodeRange(1) bodeRange(end)])
%                 title('Straight data')
%                 legend('Raw','Filtered','location','best')
%         subplot(428)
%             simpleBodephase(sysID.par.straight.sim.OL.raw ,'Hz',lw,bodeRange,'wrap');hold on;grid minor
%             simpleBodephase(sysID.par.straight.sim.OL.filt,'Hz',lw,bodeRange,'wrap');grid minor
%             simpleBodephase(G(1,1)                  ,'Hz',lw,'g--',bodeRange,'wrap');grid minor
%             xlim([bodeRange(1) bodeRange(end)])
%                 legend('Raw','Filtered','location','best')

figure(baseFig_sim+304);clf
    set(gcf,'position',[100 60 1400 450])
    sgtitle('Bode plots of identification using filtered Open Loop data')
        subplot(241)
            simpleBodemag(sysID.par.initSys.sim.OL.raw ,'Hz',lw,bodeRange);hold on;grid minor
            simpleBodemag(sysID.par.initSys.sim.OL.filt,'Hz',lw,bodeRange);grid minor
            simpleBodemag(G(1,1)                 ,'Hz',lw,'g--',bodeRange);grid minor
            xlim([bodeRange(1) bodeRange(end)])
                title('Initial system')
                legend('Using raw data','Using filtered data','location','best')
        subplot(245)
            simpleBodephase(sysID.par.initSys.sim.OL.raw ,'Hz',lw,bodeRange,'wrap');hold on;grid minor
            simpleBodephase(sysID.par.initSys.sim.OL.filt,'Hz',lw,bodeRange,'wrap');grid minor
            simpleBodephase(G(1,1)                 ,'Hz',lw,'g--',bodeRange,'wrap');grid minor
            xlim([bodeRange(1) bodeRange(end)])
                legend('Using raw data','Using filtered data','location','best')
        subplot(242)
            simpleBodemag(sysID.par.fixedOrder.sim.OL.raw ,'Hz',lw,bodeRange);hold on;grid minor
            simpleBodemag(sysID.par.fixedOrder.sim.OL.filt,'Hz',lw,bodeRange);grid minor
            simpleBodemag(G(1,1)                    ,'Hz',lw,'g--',bodeRange);grid minor
            xlim([bodeRange(1) bodeRange(end)])
                title(['Fixed order, nx=',num2str(nx)])
                legend('Using raw data','Using filtered data','location','best')
        subplot(246)
            simpleBodephase(sysID.par.fixedOrder.sim.OL.raw ,'Hz',lw,bodeRange,'wrap');hold on;grid minor
            simpleBodephase(sysID.par.fixedOrder.sim.OL.filt,'Hz',lw,bodeRange,'wrap');grid minor
            simpleBodephase(G(1,1)                    ,'Hz',lw,'g--',bodeRange,'wrap');grid minor
            xlim([bodeRange(1) bodeRange(end)])
                legend('Using raw data','Using filtered data','location','best')
        subplot(243)
            simpleBodemag(sysID.par.firstAprrox.OL.raw ,'Hz',lw,bodeRange);hold on;grid minor
            simpleBodemag(sysID.par.firstAprrox.OL.filt,'Hz',lw,bodeRange);grid minor
            simpleBodemag(G(1,1)                 ,'Hz',lw,'g--',bodeRange);grid minor
            xlim([bodeRange(1) bodeRange(end)])
                title('First order approximation')
                legend('Using raw data','Using filtered data','location','best')
        subplot(247)
            simpleBodephase(sysID.par.firstAprrox.OL.raw ,'Hz',lw,bodeRange,'wrap');hold on;grid minor
            simpleBodephase(sysID.par.firstAprrox.OL.filt,'Hz',lw,bodeRange,'wrap');grid minor
            simpleBodephase(G(1,1)                 ,'Hz',lw,'g--',bodeRange,'wrap');grid minor
            xlim([bodeRange(1) bodeRange(end)])
                legend('Using raw data','Using filtered data','location','best')
        subplot(244)
            simpleBodemag(sysID.par.straight.sim.OL.raw ,'Hz',lw,bodeRange);hold on;grid minor
            simpleBodemag(sysID.par.straight.sim.OL.filt,'Hz',lw,bodeRange);grid minor
            simpleBodemag(G(1,1)                  ,'Hz',lw,'g--',bodeRange);grid minor
            xlim([bodeRange(1) bodeRange(end)])
                title('Straight data')
                legend('Using raw data','Using filtered data','location','best')
        subplot(248)
            simpleBodephase(sysID.par.straight.sim.OL.raw ,'Hz',lw,bodeRange,'wrap');hold on;grid minor
            simpleBodephase(sysID.par.straight.sim.OL.filt,'Hz',lw,bodeRange,'wrap');grid minor
            simpleBodephase(G(1,1)                  ,'Hz',lw,'g--',bodeRange,'wrap');grid minor
            xlim([bodeRange(1) bodeRange(end)])
                legend('Using raw data','Using filtered data','location','best')
    print('./_Figures/OL/parTime/compBodeAll_1','-depsc')

%% Validate and compare the identifications in time domain, using x-step ahead prediction and a residual test

% Choose which dataset to use for the validation
    dataSet = sysID.data.OL.train.id;
    % dataSet = sysID.data.OL.val.id;
    % dataSet = sysID.data.OL.trans.id;
    
    dataSet_filt = sysID.dataFilt.OL.train.id;
    % dataSet_filt = sysID.dataFilt.OL.val.id;
    % dataSet_filt = sysID.dataFilt.OL.trans.id;

% Choose the amount of steps to use for the prediction (inf=simulation)
    % xStep = 1;
    xStep = inf;

% Run the predictions and plot them
    % test = initializedVal_OL(sysID.par.straight.sim.OL.raw,sysID.data.OL.train.id);
    OLsim_sysID_predTests_V3

% Run the residual tests and plot them
    % initializedVal_OL(sysID.par.straight.sim.OL.raw,sysID.data.OL.train.id)
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


%% Save the identification data

if saveID
    save(['../_IDdata/sysID_OL_sim_',datestr(now, 'yymmdd_HHMMSS'),'.mat'],'sysID','G','lw','tMeas','nx')
end




























