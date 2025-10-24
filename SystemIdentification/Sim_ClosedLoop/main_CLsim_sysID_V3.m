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
removeTrans = true; % [true/false] Does the transient need to be removed

% Save the identification?
saveID = true; % [true/false]

%% Make the system used for data generation
disp('Constructing system model')
main_lumpedSystem_water_V3

G = sysTMC(32,[3 6]);

%% Define some simulation variables
% fs = 1;     % [Hz]    Sampling frequency
fs = 0.1;     % [Hz]    Sampling frequency
Ts = 1/fs;  % [s]     Sampling time
% tMeas = 8; % [hours] Measurement time
tMeas = 50; % [hours] Measurement time

T0 = 23; % [degC] Inital Temperature
Tamb = 23; % [degC] Ambient Temperature

T_ref = 45; % [degC] Reference Temperature

N = tMeas*3600*fs; % [-] Measurement samples
tVec = linspace(fs,tMeas*3600,N); % [s] Time vector


%% Define the noise on the signal
v     = 0.15.*randn(N,1);
v_val = 0.15.*randn(N,1);

%% Define the identification signal
maxAmplitude = 24; % [W] Maximum exitation signal
positiveOnly = 1;  % [-] Define whether the exitation signal can be negative

disp('Generating multisine for identification')
dist     = genMultisine(fs, N, 1, maxAmplitude,     positiveOnly);
dist_val = dist;

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
    x0_CL = [T0.*ones(size(G.A,1),1); zeros(length(eig(K.tot)),1)];

% The input vector for the CL system consists of the reference temperature,
% the disturbance (/identification) signal and the ambient temperature
    u     = [T_ref.*ones(1,length(dist))    ;dist    ;ambVec];
    u_val = [T_ref.*ones(1,length(dist_val));dist_val;ambVec];

% Run the closed loop simulation, again with the previously defined noise signal added
    y_CL     = lsim(P,u    ,tVec,x0_CL)+v;
    y_CL_val = lsim(P,u_val,tVec,x0_CL)+v_val;


%% Plot the simulation results

figure(baseFig_sim+101);clf
    set(gcf,'position',[700 100 700 300])
    plot(tVec,y_CL(:,1),LineWidth=lw);grid minor
        xlabel('Time [s]')
        ylabel('Temperature [degC]')
        title(['Thermal mass temperature, open loop identification, ',num2str(tMeas),' hours'])

%% Define the in and output vector sizes

% Define the range of the data that is used for identification and validation
if removeTrans
    idx1_CL = 2.5*3600*fs;
else
    idx1_CL = 1;
end
if makeValSet
    idx2_CL = N-(N-2.5*3600*fs)*0.2;
else
    idx2_CL = N;
end

idxRange_CL_trns = 1:idx1_CL+1;
idxRange_CL      = idx1_CL:idx2_CL;
idxRange_CL_val  = (idx2_CL):N;

%% Define the in and output vectors, for easy and consistent use
% Open loop
sysID.data.CL.trans.out = y_CL(idxRange_CL_trns,1)-ambVec(idxRange_CL_trns)';
% sysID.data.CL.trans.in  = dist(idxRange_CL_trns);
sysID.data.CL.trans.in  = y_CL(idxRange_CL_trns,2)';

sysID.data.CL.train.out = y_CL(idxRange_CL,1)-ambVec(idxRange_CL)';
% sysID.data.CL.train.in  = dist(idxRange_CL);
sysID.data.CL.train.in  = y_CL(idxRange_CL,2)';

sysID.data.CL.val.out   = y_CL_val(idxRange_CL_val,1)-ambVec(idxRange_CL_val)';
% sysID.data.CL.val.in    = dist_val(idxRange_CL_val);
sysID.data.CL.val.in    = y_CL_val(idxRange_CL_val,2)';

sysID.data.CL.full.out  = y_CL(:,1)-ambVec';
% sysID.data.CL.full.in   = dist;
sysID.data.CL.full.in   = y_CL(:,2)';

% Define the time vectors
sysID.data.CL.trans.tVec = tVec(idxRange_CL_trns);
sysID.data.CL.train.tVec = tVec(idxRange_CL);
sysID.data.CL.val.tVec   = tVec(idxRange_CL_val);
sysID.data.CL.full.tVec  = tVec;

% tVec_CL      = tVec(idxRange_CL);

%% Filter the data vectors to remove some of the effect of noise
% High order lowpass filter, to filter away all the noise above the 
% multisine exitation frequency
    LPfilt = orderLP(0.01,2,0,0,0.7)^2;
    % [LPfiltB,LPfiltA] = butter(4,0.005);
    % LPfilt = tf(LPfiltB,LPfiltA,1);

% Actually filter the data
    % Open loop data
    sysID.dataFilt.CL.trans.out = lsim(LPfilt,sysID.data.CL.trans.out-sysID.data.CL.trans.out(1),sysID.data.CL.trans.tVec)'+sysID.data.CL.trans.out(1);
    sysID.dataFilt.CL.trans.in = lsim(LPfilt,sysID.data.CL.trans.in-sysID.data.CL.trans.in(1),sysID.data.CL.trans.tVec)'+sysID.data.CL.trans.in(1);

    sysID.dataFilt.CL.train.out = lsim(LPfilt,sysID.data.CL.train.out-sysID.data.CL.train.out(1),sysID.data.CL.train.tVec)'+sysID.data.CL.train.out(1);
    sysID.dataFilt.CL.train.in = lsim(LPfilt,sysID.data.CL.train.in-sysID.data.CL.train.in(1),sysID.data.CL.train.tVec)'+sysID.data.CL.train.in(1);

    sysID.dataFilt.CL.val.out = lsim(LPfilt,sysID.data.CL.val.out-sysID.data.CL.val.out(1),sysID.data.CL.val.tVec)'+sysID.data.CL.val.out(1);
    sysID.dataFilt.CL.val.in = lsim(LPfilt,sysID.data.CL.val.in-sysID.data.CL.val.in(1),sysID.data.CL.val.tVec)'+sysID.data.CL.val.in(1);

    sysID.dataFilt.CL.full.out = lsim(LPfilt,sysID.data.CL.full.out-sysID.data.CL.full.out(1),sysID.data.CL.full.tVec)'+sysID.data.CL.full.out(1);
    sysID.dataFilt.CL.full.in = lsim(LPfilt,sysID.data.CL.full.in-sysID.data.CL.full.in(1),sysID.data.CL.full.tVec)'+sysID.data.CL.full.in(1);
    
%% FRF measurements using the new state vector
% Make the FRF's using the raw data
    [sysID.nonPar.trd.CL.raw, ~] = makeOpenLoopFRF_sysIdent(sysID.data.CL.train.out, sysID.data.CL.train.in(1,:), fs);

% Make the FRF's using the filtered data
    [sysID.nonPar.trd.CL.filt, ~] = makeOpenLoopFRF_sysIdent(sysID.dataFilt.CL.train.out, sysID.dataFilt.CL.train.in(1,:), fs);

% Define and do the non-parametric identification using the LPM, using both
% OL and CL data
    polyOrder = 6; % The order of the polynomial that is fitted to the data
    locality  = 8; % Amount of points (pos & neg) to consider around apprx freq
    
        disp('Starting CL LPM identification')
    [sysID.nonPar.lpm.CL.raw,~]  = sysID_LPM(sysID.data.CL.train.out     ,sysID.data.CL.train.in(1,:)     ,fs,locality,polyOrder);
    [sysID.nonPar.lpm.CL.filt,~] = sysID_LPM(sysID.dataFilt.CL.train.out,sysID.dataFilt.CL.train.in(1,:),fs,locality,polyOrder);
        disp('Finished CL LPM identification')
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
            simpleBodemag(sysID.nonPar.trd.CL.raw,'Hz',lw);hold on;grid minor
            simpleBodemag(sysID.nonPar.lpm.CL.raw,'Hz',lw);grid minor
            simpleBodemag(G(1,1),'Hz',lw,'g--');grid minor
            xlim([1e-4 0.01])
                title(['Using unfiltered data'])
                legend('Traditional','LPM','Model','location','best')
        subplot(223)
            simpleBodephase(sysID.nonPar.trd.CL.raw,'Hz',lw,'wrap');hold on;grid minor
            simpleBodephase(sysID.nonPar.lpm.CL.raw,'Hz',lw,'wrap');grid minor
            simpleBodephase(G(1,1),'Hz',lw,'g--','wrap');grid minor
            xlim([1e-4 0.01])
                legend('Traditional','LPM','Model','location','best')
        subplot(222)
            simpleBodemag(sysID.nonPar.trd.CL.filt,'Hz',lw);hold on;grid minor
            simpleBodemag(sysID.nonPar.lpm.CL.filt,'Hz',lw);grid minor
            simpleBodemag(G(1,1),'Hz',lw,'g--');grid minor
            xlim([1e-4 0.01])
                title(['Using filtered data'])
                legend('Traditional','LPM','Model','location','best')
        subplot(224)
            simpleBodephase(sysID.nonPar.trd.CL.filt,'Hz',lw,'wrap');hold on;grid minor
            simpleBodephase(sysID.nonPar.lpm.CL.filt,'Hz',lw,'wrap');grid minor
            simpleBodephase(G(1,1),'Hz',lw,'g--','wrap');grid minor
            xlim([1e-4 0.01])
                legend('Traditional','LPM','Model','location','best')

%% First order parametric approximation, using time data
    sysID.par.firstAprrox.CL.raw  = step_sysID(sysID.data.CL.full.in',zeros(size(sysID.data.CL.full.in))',sysID.data.CL.full.out,sysID.data.CL.full.tVec,10000*fs);
    sysID.par.firstAprrox.CL.filt = step_sysID(sysID.dataFilt.CL.full.in',zeros(size(sysID.dataFilt.CL.full.in))',sysID.dataFilt.CL.full.out',sysID.data.CL.full.tVec,10000*fs);

%% Visualization of the first order approximation
bodeRange = logspace(-6,-1,100);

figure(baseFig_sim+203);clf
    set(gcf,'position',[500 100 900 500])
    sgtitle(['Heater to thermal mass temperature, first order approximation, ',num2str(tMeas),' hours'])
        subplot(211)
            simpleBodemag(sysID.par.firstAprrox.CL.raw,'Hz',lw ,bodeRange);hold on;grid minor
            simpleBodemag(sysID.par.firstAprrox.CL.filt,'Hz',lw,bodeRange,'-.');grid minor
            simpleBodemag(G(1,1),'Hz',lw,'g--',bodeRange);grid minor
                legend('Gotten from unfiltered data','Gotten from filtered data','Model','location','best')
        subplot(212)
            simpleBodephase(sysID.par.firstAprrox.CL.raw,'Hz',lw,'wrap' ,bodeRange);hold on;grid minor
            simpleBodephase(sysID.par.firstAprrox.CL.filt,'Hz',lw,'wrap',bodeRange,'-.');grid minor
            simpleBodephase(G(1,1),'Hz',lw,'g--','wrap',bodeRange);grid minor
                legend('Gotten from unfiltered data','Gotten from filtered data','Model','location','best')

%% Parametric system identification, using time data, pre-requisites

% Definitions for identification
    nx = 4; % Model order for fixed order identification

% nx order initial system
    initA = sysID.par.firstAprrox.CL.raw.A.*eye(nx);
    initB = zeros(nx,1);initB(1) = sysID.par.firstAprrox.CL.raw.B;
    initC = zeros(1,nx);initC(1) = sysID.par.firstAprrox.CL.raw.C;
    
    init_sys = idss(initA,initB,initC,0,zeros(nx,1),zeros(nx,1),0);
    init_sys.Structure.B.Free = ones(nx,1);
    init_sys.Structure.C.Free = ones(1,nx);
    init_sys.Structure.K.Free = zeros(nx,1);

% Define the data as iddata's
    sysID.data.CL.train.id      = iddata(sysID.data.CL.train.out          ,sysID.data.CL.train.in(:,:)'    ,Ts);
    sysID.dataFilt.CL.train.id = iddata(sysID.dataFilt.CL.train.out'    ,sysID.dataFilt.CL.train.in'    ,Ts);

% Define the data that is not used in the identification as validation data

if makeValSet
    sysID.data.CL.val.id      = iddata(sysID.data.CL.val.out      ,sysID.data.CL.val.in',Ts);
    sysID.dataFilt.CL.val.id = iddata(sysID.dataFilt.CL.val.out',sysID.data.CL.val.in',Ts);
else
    sysID.data.CL.val.id      = iddata(1,1,Ts);
    sysID.dataFilt.CL.val.id = iddata(1,1,Ts);
end
if removeTrans
    sysID.data.CL.trans.id      = iddata(sysID.data.CL.trans.out      ,sysID.data.CL.trans.in',Ts);
    sysID.dataFilt.CL.trans.id = iddata(sysID.dataFilt.CL.trans.out',sysID.data.CL.trans.in',Ts);
else
    sysID.data.CL.trans.id      = iddata(1,1,Ts);
    sysID.dataFilt.CL.trans.id = iddata(1,1,Ts);
end

optSS_sim = ssestOptions('Focus','Simulation');
optSS_prd = ssestOptions('Focus','Prediction');

%% Parametric system identification, using time data

% All tested parametric identification options using OL data, using prediction focus
        disp(['Open Loop identification using a fixed order (nx=',num2str(nx),'): in progress'])
    sysID.par.fixedOrder.sim.CL.raw  = ssest(sysID.data.CL.train.id,nx,optSS_sim);
    sysID.par.fixedOrder.sim.CL.filt = ssest(sysID.dataFilt.CL.train.id,nx,optSS_sim);
    sysID.par.firstAprrox.CL.raw     = ssest(sysID.data.CL.train.id,1,optSS_sim);
    sysID.par.firstAprrox.CL.filt    = ssest(sysID.dataFilt.CL.train.id,1,optSS_sim);
        disp(['Open Loop identification using a fixed order (nx=',num2str(nx),'): done'])

        disp('Open Loop identification using an initial system: in progress')
    sysID.par.initSys.sim.CL.raw  = ssest(sysID.data.CL.train.id    ,sysID.par.firstAprrox.CL.raw,optSS_sim);
    sysID.par.initSys.sim.CL.filt = ssest(sysID.dataFilt.CL.train.id,sysID.par.firstAprrox.CL.filt,optSS_sim);
    % sysID.par.initSys.sim.CL.raw  = ssest(sysID.data.CL.train.id    ,sysID.par.fixedOrder.prd.CL.raw,optSS_sim);
    % sysID.par.initSys.sim.CL.filt = ssest(sysID.dataFilt.CL.train.id,sysID.par.fixedOrder.prd.CL.raw,optSS_sim);
        disp('Open Loop identification using an initial system: done')

% All tested parametric identification options using CL data, using prediction focus
        disp(['Open Loop identification using a fixed order (nx=',num2str(nx),'): in progress'])
    sysID.par.fixedOrder.prd.CL.raw  = ssest(sysID.data.CL.train.id,nx,optSS_prd);
    sysID.par.fixedOrder.prd.CL.filt = ssest(sysID.dataFilt.CL.train.id,nx,optSS_prd);
    sysID.par.firstAprrox.CL.raw     = ssest(sysID.data.CL.train.id,1,optSS_prd);
    sysID.par.firstAprrox.CL.filt    = ssest(sysID.dataFilt.CL.train.id,1,optSS_prd);
        disp(['Open Loop identification using a fixed order (nx=',num2str(nx),'): done'])

        disp('Open Loop identification using an initial system: in progress')
    sysID.par.initSys.prd.CL.raw  = ssest(sysID.data.CL.train.id    ,sysID.par.firstAprrox.CL.raw,optSS_prd);
    sysID.par.initSys.prd.CL.filt = ssest(sysID.dataFilt.CL.train.id,sysID.par.firstAprrox.CL.filt,optSS_prd);
    % sysID.par.initSys.prd.CL.raw  = ssest(sysID.data.CL.train.id    ,sysID.par.fixedOrder.prd.CL.raw,optSS_prd);
    % sysID.par.initSys.prd.CL.filt = ssest(sysID.dataFilt.CL.train.id,sysID.par.fixedOrder.prd.CL.raw,optSS_prd);
        disp('Open Loop identification using an initial system: done')

        
        disp('Open Loop identification using the raw data: in progress')
    sysID.par.straight.sim.CL.raw  = ssest(sysID.data.CL.train.id);
    sysID.par.straight.sim.CL.filt = ssest(sysID.dataFilt.CL.train.id);
        disp('Open Loop identification using the raw data: done')

%% Visualize and compare the identifications in frequency domain
bodeRange = logspace(-6,-1,1000);

figure(baseFig_sim+301);clf
    % set(gcf,'position',[500 100 700 450])
    sgtitle('Bode plots of identification using raw Open Loop data')
        subplot(211)
            simpleBodemag(sysID.par.initSys.sim.CL.raw   ,'Hz',lw,bodeRange);hold on;grid minor
            simpleBodemag(sysID.par.fixedOrder.sim.CL.raw,'Hz',lw,bodeRange);grid minor
            simpleBodemag(sysID.par.firstAprrox.CL.raw   ,'Hz',lw,bodeRange);grid minor
            simpleBodemag(sysID.par.straight.sim.CL.raw  ,'Hz',lw,bodeRange);grid minor
            simpleBodemag(G(1,1)                   ,'Hz',lw,'g--',bodeRange);grid minor
                legend('InitSysIdent','Fixed order','First order approx','Straight data','Model','location','best')
        subplot(212)
            simpleBodephase(sysID.par.initSys.sim.CL.raw   ,'Hz',lw,bodeRange,'wrap');hold on;grid minor
            simpleBodephase(sysID.par.fixedOrder.sim.CL.raw,'Hz',lw,bodeRange,'wrap');grid minor
            simpleBodephase(sysID.par.firstAprrox.CL.raw   ,'Hz',lw,bodeRange,'wrap');grid minor
            simpleBodephase(sysID.par.straight.sim.CL.raw  ,'Hz',lw,bodeRange,'wrap');grid minor
            simpleBodephase(G(1,1)                   ,'Hz',lw,'g--',bodeRange,'wrap');grid minor
                legend('InitSysIdent','Fixed order','First order approx','Straight data','Model','location','best')

figure(baseFig_sim+302);clf
    % set(gcf,'position',[500 100 700 450])
    sgtitle('Bode plots of identification using raw Open Loop data')
        subplot(211)
            simpleBodemag(sysID.par.initSys.prd.CL.raw   ,'Hz',lw,bodeRange);hold on;grid minor
            simpleBodemag(sysID.par.fixedOrder.prd.CL.raw,'Hz',lw,bodeRange);grid minor
            simpleBodemag(sysID.par.firstAprrox.CL.raw   ,'Hz',lw,bodeRange);grid minor
            simpleBodemag(sysID.par.straight.sim.CL.raw  ,'Hz',lw,bodeRange);grid minor
            simpleBodemag(G(1,1)                   ,'Hz',lw,'g--',bodeRange);grid minor
                legend('InitSysIdent','Fixed order','First order approx','Straight data','Model','location','best')
        subplot(212)
            simpleBodephase(sysID.par.initSys.prd.CL.raw   ,'Hz',lw,bodeRange,'wrap');hold on;grid minor
            simpleBodephase(sysID.par.fixedOrder.prd.CL.raw,'Hz',lw,bodeRange,'wrap');grid minor
            simpleBodephase(sysID.par.firstAprrox.CL.raw   ,'Hz',lw,bodeRange,'wrap');grid minor
            simpleBodephase(sysID.par.straight.sim.CL.raw  ,'Hz',lw,bodeRange,'wrap');grid minor
            simpleBodephase(G(1,1)                   ,'Hz',lw,'g--',bodeRange,'wrap');grid minor
                legend('InitSysIdent','Fixed order','First order approx','Straight data','Model','location','best')

figure(baseFig_sim+303);clf
    % set(gcf,'position',[500 100 700 450])
    sgtitle('Bode plots of identification using filtered Open Loop data')
        subplot(211)
            simpleBodemag(sysID.par.initSys.sim.CL.filt   ,'Hz',lw,bodeRange);hold on;grid minor
            simpleBodemag(sysID.par.fixedOrder.sim.CL.filt,'Hz',lw,bodeRange);grid minor
            simpleBodemag(sysID.par.firstAprrox.CL.filt   ,'Hz',lw,bodeRange);grid minor
            simpleBodemag(sysID.par.straight.sim.CL.filt  ,'Hz',lw,bodeRange);grid minor
            simpleBodemag(G(1,1)                    ,'Hz',lw,'g--',bodeRange);grid minor
                legend('InitSysIdent','Fixed order','First order approx','Straight data','Model','location','best')
        subplot(212)
            simpleBodephase(sysID.par.initSys.sim.CL.filt   ,'Hz',lw,bodeRange,'wrap');hold on;grid minor
            simpleBodephase(sysID.par.fixedOrder.sim.CL.filt,'Hz',lw,bodeRange,'wrap');grid minor
            simpleBodephase(sysID.par.firstAprrox.CL.filt   ,'Hz',lw,bodeRange,'wrap');grid minor
            simpleBodephase(sysID.par.straight.sim.CL.filt  ,'Hz',lw,bodeRange,'wrap');grid minor
            simpleBodephase(G(1,1)                    ,'Hz',lw,'g--',bodeRange,'wrap');grid minor
                legend('InitSysIdent','Fixed order','First order approx','Straight data','Model','location','best')

% figure(baseFig_sim+303);clf
%     set(gcf,'position',[100 60 1200 800])
%     sgtitle('Bode plots of identification using filtered Open Loop data')
%         subplot(421)
%             simpleBodemag(sysID.par.initSys.sim.CL.raw ,'Hz',lw,bodeRange);hold on;grid minor
%             simpleBodemag(sysID.par.initSys.sim.CL.filt,'Hz',lw,bodeRange);grid minor
%             simpleBodemag(G(1,1)                 ,'Hz',lw,'g--',bodeRange);grid minor
%             xlim([bodeRange(1) bodeRange(end)])
%                 title('Initial system')
%                 legend('Raw','Filtered','location','best')
%         subplot(423)
%             simpleBodephase(sysID.par.initSys.sim.CL.raw ,'Hz',lw,bodeRange,'wrap');hold on;grid minor
%             simpleBodephase(sysID.par.initSys.sim.CL.filt,'Hz',lw,bodeRange,'wrap');grid minor
%             simpleBodephase(G(1,1)                 ,'Hz',lw,'g--',bodeRange,'wrap');grid minor
%             xlim([bodeRange(1) bodeRange(end)])
%                 legend('Raw','Filtered','location','best')
%         subplot(422)
%             simpleBodemag(sysID.par.fixedOrder.sim.CL.raw ,'Hz',lw,bodeRange);hold on;grid minor
%             simpleBodemag(sysID.par.fixedOrder.sim.CL.filt,'Hz',lw,bodeRange);grid minor
%             simpleBodemag(G(1,1)                    ,'Hz',lw,'g--',bodeRange);grid minor
%             xlim([bodeRange(1) bodeRange(end)])
%                 title(['Fixed order, nx=',num2str(nx)])
%                 legend('Raw','Filtered','location','best')
%         subplot(424)
%             simpleBodephase(sysID.par.fixedOrder.sim.CL.raw ,'Hz',lw,bodeRange,'wrap');hold on;grid minor
%             simpleBodephase(sysID.par.fixedOrder.sim.CL.filt,'Hz',lw,bodeRange,'wrap');grid minor
%             simpleBodephase(G(1,1)                    ,'Hz',lw,'g--',bodeRange,'wrap');grid minor
%             xlim([bodeRange(1) bodeRange(end)])
%                 legend('Raw','Filtered','location','best')
%         subplot(425)
%             simpleBodemag(sysID.par.firstAprrox.CL.raw ,'Hz',lw,bodeRange);hold on;grid minor
%             simpleBodemag(sysID.par.firstAprrox.CL.filt,'Hz',lw,bodeRange);grid minor
%             simpleBodemag(G(1,1)                 ,'Hz',lw,'g--',bodeRange);grid minor
%             xlim([bodeRange(1) bodeRange(end)])
%                 title('First order approximation')
%                 legend('Raw','Filtered','location','best')
%         subplot(427)
%             simpleBodephase(sysID.par.firstAprrox.CL.raw ,'Hz',lw,bodeRange,'wrap');hold on;grid minor
%             simpleBodephase(sysID.par.firstAprrox.CL.filt,'Hz',lw,bodeRange,'wrap');grid minor
%             simpleBodephase(G(1,1)                 ,'Hz',lw,'g--',bodeRange,'wrap');grid minor
%             xlim([bodeRange(1) bodeRange(end)])
%                 legend('Raw','Filtered','location','best')
%         subplot(426)
%             simpleBodemag(sysID.par.straight.sim.CL.raw ,'Hz',lw,bodeRange);hold on;grid minor
%             simpleBodemag(sysID.par.straight.sim.CL.filt,'Hz',lw,bodeRange);grid minor
%             simpleBodemag(G(1,1)                  ,'Hz',lw,'g--',bodeRange);grid minor
%             xlim([bodeRange(1) bodeRange(end)])
%                 title('Straight data')
%                 legend('Raw','Filtered','location','best')
%         subplot(428)
%             simpleBodephase(sysID.par.straight.sim.CL.raw ,'Hz',lw,bodeRange,'wrap');hold on;grid minor
%             simpleBodephase(sysID.par.straight.sim.CL.filt,'Hz',lw,bodeRange,'wrap');grid minor
%             simpleBodephase(G(1,1)                  ,'Hz',lw,'g--',bodeRange,'wrap');grid minor
%             xlim([bodeRange(1) bodeRange(end)])
%                 legend('Raw','Filtered','location','best')

figure(baseFig_sim+304);clf
    set(gcf,'position',[100 60 1400 450])
    sgtitle('Bode plots of identification using filtered Open Loop data')
        subplot(241)
            simpleBodemag(sysID.par.initSys.sim.CL.raw ,'Hz',lw,bodeRange);hold on;grid minor
            simpleBodemag(sysID.par.initSys.sim.CL.filt,'Hz',lw,bodeRange);grid minor
            simpleBodemag(G(1,1)                 ,'Hz',lw,'g--',bodeRange);grid minor
            xlim([bodeRange(1) bodeRange(end)])
                title('Initial system')
                legend('Using raw data','Using filtered data','location','best')
        subplot(245)
            simpleBodephase(sysID.par.initSys.sim.CL.raw ,'Hz',lw,bodeRange,'wrap');hold on;grid minor
            simpleBodephase(sysID.par.initSys.sim.CL.filt,'Hz',lw,bodeRange,'wrap');grid minor
            simpleBodephase(G(1,1)                 ,'Hz',lw,'g--',bodeRange,'wrap');grid minor
            xlim([bodeRange(1) bodeRange(end)])
                legend('Using raw data','Using filtered data','location','best')
        subplot(242)
            simpleBodemag(sysID.par.fixedOrder.sim.CL.raw ,'Hz',lw,bodeRange);hold on;grid minor
            simpleBodemag(sysID.par.fixedOrder.sim.CL.filt,'Hz',lw,bodeRange);grid minor
            simpleBodemag(G(1,1)                    ,'Hz',lw,'g--',bodeRange);grid minor
            xlim([bodeRange(1) bodeRange(end)])
                title(['Fixed order, nx=',num2str(nx)])
                legend('Using raw data','Using filtered data','location','best')
        subplot(246)
            simpleBodephase(sysID.par.fixedOrder.sim.CL.raw ,'Hz',lw,bodeRange,'wrap');hold on;grid minor
            simpleBodephase(sysID.par.fixedOrder.sim.CL.filt,'Hz',lw,bodeRange,'wrap');grid minor
            simpleBodephase(G(1,1)                    ,'Hz',lw,'g--',bodeRange,'wrap');grid minor
            xlim([bodeRange(1) bodeRange(end)])
                legend('Using raw data','Using filtered data','location','best')
        subplot(243)
            simpleBodemag(sysID.par.firstAprrox.CL.raw ,'Hz',lw,bodeRange);hold on;grid minor
            simpleBodemag(sysID.par.firstAprrox.CL.filt,'Hz',lw,bodeRange);grid minor
            simpleBodemag(G(1,1)                 ,'Hz',lw,'g--',bodeRange);grid minor
            xlim([bodeRange(1) bodeRange(end)])
                title('First order approximation')
                legend('Using raw data','Using filtered data','location','best')
        subplot(247)
            simpleBodephase(sysID.par.firstAprrox.CL.raw ,'Hz',lw,bodeRange,'wrap');hold on;grid minor
            simpleBodephase(sysID.par.firstAprrox.CL.filt,'Hz',lw,bodeRange,'wrap');grid minor
            simpleBodephase(G(1,1)                 ,'Hz',lw,'g--',bodeRange,'wrap');grid minor
            xlim([bodeRange(1) bodeRange(end)])
                legend('Using raw data','Using filtered data','location','best')
        subplot(244)
            simpleBodemag(sysID.par.straight.sim.CL.raw ,'Hz',lw,bodeRange);hold on;grid minor
            simpleBodemag(sysID.par.straight.sim.CL.filt,'Hz',lw,bodeRange);grid minor
            simpleBodemag(G(1,1)                  ,'Hz',lw,'g--',bodeRange);grid minor
            xlim([bodeRange(1) bodeRange(end)])
                title('Straight data')
                legend('Using raw data','Using filtered data','location','best')
        subplot(248)
            simpleBodephase(sysID.par.straight.sim.CL.raw ,'Hz',lw,bodeRange,'wrap');hold on;grid minor
            simpleBodephase(sysID.par.straight.sim.CL.filt,'Hz',lw,bodeRange,'wrap');grid minor
            simpleBodephase(G(1,1)                  ,'Hz',lw,'g--',bodeRange,'wrap');grid minor
            xlim([bodeRange(1) bodeRange(end)])
                legend('Using raw data','Using filtered data','location','best')

%% Validate and compare the identifications in time domain, using x-step ahead prediction and a residual test

% Choose which dataset to use for the validation
    % dataSet = sysID.data.CL.train.id;
    dataSet = sysID.data.CL.val.id;
    % dataSet = sysID.data.CL.trans.id;
    
    % dataSet_filt = sysID.dataFilt.CL.train.id;
    dataSet_filt = sysID.dataFilt.CL.val.id;
    % dataSet_filt = sysID.dataFilt.CL.trans.id;

% Choose the amount of steps to use for the prediction (inf=simulation)
    % xStep = 1;
    xStep = inf;

% Run the predictions and plot them
    CLsim_sysID_predTests_V3

% Run the residual tests and plot them
    % OLsim_sysID_resTests_V3


%% Frequency domain identification

nx_freq = 4;

CL_freqDat_trd      = idfrd(squeeze(sysID.nonPar.trd.CL.raw.ResponseData) ,sysID.nonPar.trd.CL.raw.Frequency ,1);
CL_freqDat_lpm      = idfrd(squeeze(sysID.nonPar.lpm.CL.raw.ResponseData) ,sysID.nonPar.lpm.CL.raw.Frequency ,1);
CL_freqDat_trd_filt = idfrd(squeeze(sysID.nonPar.trd.CL.filt.ResponseData),sysID.nonPar.trd.CL.filt.Frequency,1);
CL_freqDat_lpm_filt = idfrd(squeeze(sysID.nonPar.lpm.CL.filt.ResponseData),sysID.nonPar.lpm.CL.filt.Frequency,1);

sysID.par.freqTRD.sim.CL.raw = ssest(CL_freqDat_trd,nx_freq,optSS_sim);
sysID.par.freqLPM.sim.CL.raw = ssest(CL_freqDat_lpm,nx_freq,optSS_sim);
% sysID.par.freqTRD.sim.CL.raw  = ssest(CL_freqDat_trd     ,init_sys,optSS_sim);
% sysID.par.freqLPM.sim.CL.raw  = ssest(CL_freqDat_lpm     ,init_sys,optSS_sim);
% sysID.par.freqTRD.sim.CL.raw  = ssest(CL_freqDat_trd     ,sysID.par.fixedOrder.sim.CL.raw,optSS_sim);
% sysID.par.freqLPM.sim.CL.raw  = ssest(CL_freqDat_lpm     ,sysID.par.fixedOrder.sim.CL.raw,optSS_sim);
% sysID.par.freqTRD.sim.CL.filt = ssest(CL_freqDat_trd_filt,sysID.par.fixedOrder.sim.CL.raw,optSS_sim);
% sysID.par.freqLPM.sim.CL.filt = ssest(CL_freqDat_lpm_filt,sysID.par.fixedOrder.sim.CL.raw,optSS_sim);

%%
figure(baseFig_sim+701);clf
    % bode(sysID.par.freqTRD.sim.CL.raw,sysID.par.freqLPM.sim.CL.raw,sysID.par.freqTRD.sim.CL.filt,sysID.par.freqLPM.sim.CL.filt,G(1,1),'g--',opt);grid minor
    bode(sysID.par.freqTRD.sim.CL.raw,sysID.par.freqLPM.sim.CL.raw,G(1,1),'g--',opt);grid minor


%% Time based validation using different (CL) data set
    
CLrange = 15000:length(datCL.Watt);
% CLrange = 1:length(datCL.Watt);

% datCL.Watt = datCL.Watt.*5;


err_initSys     = lsim(sysID.par.initSys.sim.CL.filt   ,datCL.Watt,datCL.tVec)+datCL.tempAmb'-datCL.tempTM';
err_fixedOrder  = lsim(sysID.par.fixedOrder.sim.CL.filt,datCL.Watt,datCL.tVec)+datCL.tempAmb'-datCL.tempTM';
err_firstAprrox = lsim(sysID.par.firstAprrox.CL.filt   ,datCL.Watt,datCL.tVec)+datCL.tempAmb'-datCL.tempTM';
err_straight    = lsim(sysID.par.straight.sim.CL.filt  ,datCL.Watt,datCL.tVec)+datCL.tempAmb'-datCL.tempTM';
    
sig_initSys     = lsim(sysID.par.initSys.sim.CL.filt   ,datCL.Watt,datCL.tVec)+datCL.tempAmb';
sig_fixedOrder  = lsim(sysID.par.fixedOrder.sim.CL.filt,datCL.Watt,datCL.tVec)+datCL.tempAmb';
sig_firstAprrox = lsim(sysID.par.firstAprrox.CL.filt   ,datCL.Watt,datCL.tVec)+datCL.tempAmb';
sig_straight    = lsim(sysID.par.straight.sim.CL.filt  ,datCL.Watt,datCL.tVec)+datCL.tempAmb';

figure(baseFig_sim+651);clf
    subplot(211);hold on;grid minor
        plot(datCL.tVec(CLrange),sig_initSys(CLrange))
        plot(datCL.tVec(CLrange),sig_fixedOrder(CLrange))
        plot(datCL.tVec(CLrange),sig_firstAprrox(CLrange))
        plot(datCL.tVec(CLrange),sig_straight(CLrange))
        plot(datCL.tVec(CLrange),datCL.tempTM(CLrange))
            xlabel('Time [s]')
            ylabel('Temperature [degC]')
            title('Comparison of measured and simulated data using identified model')
            legend('Simulated data using identified model (initSys)', ...
                   'Simulated data using identified model (fixed order)', ...
                   'Simulated data using identified model (firstApprox)', ...
                   'Simulated data using identified model (Straight)', ...
                   'Measured data','location','southeast')
    subplot(212);hold on;grid minor
        plot(datCL.tVec(CLrange),detrend(err_initSys(CLrange),1))
        plot(datCL.tVec(CLrange),detrend(err_fixedOrder(CLrange),1))
        plot(datCL.tVec(CLrange),detrend(err_firstAprrox(CLrange),1))
        plot(datCL.tVec(CLrange),detrend(err_straight(CLrange),1))
            xlabel('Time [s]')
            ylabel('Temperature diff. [degC]')
            title('Difference between measured and simulated data using identified model')
            legend(['Simulated data using identified model (initSys), RMSE: ',num2str(rms(err_initSys),2),'degC'], ...
                   ['Simulated data using identified model (fixed order), RMSE: ',num2str(rms(err_fixedOrder),2),'degC'], ...
                   ['Simulated data using identified model (firstApprox), RMSE: ',num2str(rms(err_firstAprrox),2),'degC'], ...
                   ['Simulated data using identified model (Straight), RMSE: ',num2str(rms(err_straight),2),'degC'], ...
               'location','southeast')

% datCL.Watt = datCL.Watt./5;

[PSD_initsys    ,CAS_initsys    ,freqVec_initsys]     = fftCas_V2(err_initSys(CLrange)    ,datCL.Ts);
[PSD_fixedOrder ,CAS_fixedOrder ,freqVec_fixedOrder]  = fftCas_V2(err_fixedOrder(CLrange) ,datCL.Ts);
[PSD_firstAprrox,CAS_firstAprrox,freqVec_firstAprrox] = fftCas_V2(err_firstAprrox(CLrange),datCL.Ts);
[PSD_straight   ,CAS_straight   ,freqVec_straight]    = fftCas_V2(err_straight(CLrange)   ,datCL.Ts);
    
figure(baseFig_sim+652);clf
    subplot(211)
        semilogx(freqVec_initsys    ,db(squeeze(PSD_initsys)));grid minor;hold on
        semilogx(freqVec_fixedOrder ,db(squeeze(PSD_fixedOrder)))
        semilogx(freqVec_firstAprrox,db(squeeze(PSD_firstAprrox)))
        semilogx(freqVec_straight   ,db(squeeze(PSD_straight)))
            xlabel('Frequency [Hz]')
            ylabel('Power [$degC^2/Hz$]')
            title('PSD of simulation error')
            legend('initSys', ...
                   'fixed order', ...
                   'firstApprox', ...
                   'Straight','location','southeast')
    subplot(212)
        semilogx(freqVec_initsys    ,CAS_initsys    );grid minor;hold on
        semilogx(freqVec_fixedOrder ,CAS_fixedOrder )
        semilogx(freqVec_firstAprrox,CAS_firstAprrox)
        semilogx(freqVec_straight   ,CAS_straight   )
            xlabel('Frequency [Hz]')
            ylabel('RMSE [degC]')
            title('CAS of simulation error')


%% Save the identification data

% if saveID
%     save(['../_IDdata/sysID_CL_sim_',datestr(now, 'yymmdd_HHMMSS'),'.mat'],'sysID','G','lw','tMeas','nx')
% end




























