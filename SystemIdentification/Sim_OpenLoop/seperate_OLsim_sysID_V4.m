%———————————————————————————————————————————————————————————————————————————
%  ██████╗ ██████╗ ███████╗███╗  ██╗     ██╗      ██████╗  ██████╗ ██████╗ 
% ██╔═══██╗██╔══██╗██╔════╝████╗ ██║     ██║     ██╔═══██╗██╔═══██╗██╔══██╗
% ██║   ██║██████╔╝███████╗██╔██╗██║     ██║     ██║   ██║██║   ██║██████╔╝
% ██║   ██║██╔═══╝ ██╔════╝██║╚████║     ██║     ██║   ██║██║   ██║██╔═══╝ 
% ╚██████╔╝██║     ███████╗██║ ╚███║     ███████╗╚██████╔╝╚██████╔╝██║     
%  ╚═════╝ ╚═╝     ╚══════╝╚═╝  ╚══╝     ╚══════╝ ╚═════╝  ╚═════╝ ╚═╝     
%
%  ██████╗██╗███╗   ███╗██╗   ██╗██╗      █████╗ ████████╗██╗ ██████╗ ███╗  ██╗
% ██╔════╝██║████╗ ████║██║   ██║██║     ██╔══██╗╚══██╔══╝██║██╔═══██╗████╗ ██║
% ╚█████╗ ██║██╔████╔██║██║   ██║██║     ██║  ██║   ██║   ██║██║   ██║██╔██╗██║
%  ╚═══██╗██║██║╚██╔╝██║██║   ██║██║     ███████║   ██║   ██║██║   ██║██║╚████║
% ██████╔╝██║██║ ╚═╝ ██║╚██████╔╝███████╗██╔══██║   ██║   ██║╚██████╔╝██║ ╚███║
% ╚═════╝ ╚═╝╚═╝     ╚═╝ ╚═════╝ ╚══════╝╚═╝  ╚═╝   ╚═╝   ╚═╝ ╚═════╝ ╚═╝  ╚══╝
%———————————————————————————————————————————————————————————————————————————————

%% Run simulation and format data, only if new data is requested
% Always uses the pre-defined disturbance and noise signal, since these
% might be needed to stay constant over several simulations
if newData
    [sysID.OLdata, sysID.OLdataFilt] = makeOLdata(G,tMeas,dist,v,Ts);
end

%% First order parametric approximation, using direct time data
    sysID.par.OL.time.firstApprox.raw  = firstOrderApprox_OL(sysID.OLdata.train.id);
    sysID.par.OL.time.firstApprox.filt = firstOrderApprox_OL(sysID.OLdataFilt.train.id);

%% FRF measurements
if compTrans

disp('Starting OL FRF identification')
% Make the FRF's using the raw uncompensated data
    [sysID.nonPar.OL.frf.trd.raw, ~] = makeOpenLoopFRF_sysIdent(sysID.OLdata.train.y, sysID.OLdata.train.u, fs);

% Make the FRF's using the filtered uncompensated data
    [sysID.nonPar.OL.frf.trd.filt, ~] = makeOpenLoopFRF_sysIdent(sysID.OLdataFilt.train.y, sysID.OLdata.train.u, fs);
disp('Finished OL FRF identification')

% Define and do the non-parametric identification using the LPM
    polyOrder = 6; % The order of the polynomial that is fitted to the data
    locality  = 8; % Amount of points (pos & neg) to consider around apprx freq
    
    disp('Starting OL LPM identification')
        [sysID.nonPar.OL.frf.lpm.raw,~]  = sysID_LPM(sysID.OLdata.train.y     ,sysID.OLdata.train.u'     ,fs,locality,polyOrder);
        [sysID.nonPar.OL.frf.lpm.filt,~] = sysID_LPM(sysID.OLdataFilt.train.y,sysID.OLdata.train.u',fs,locality,polyOrder);
    disp('Finished OL LPM identification')

if plotAll
    OL_sim_validate_nonParam
end

% Using an exponential function fit, remove the transient from the data
    % Format the data
        trainData     = sysID.OLdata.train;
        trainDataFilt = sysID.OLdataFilt.train;
        timeData      = sysID.OLdata.train.tVec;
        
    % Calculate the steady state value of the measurement
        DCgain = mean(trainData.y(end-1000:end));
    
    % Make an exponential fit using Matlab's 'fit.m' function
        detrendOut = fit(timeData',trainData.y-DCgain,'exp1');
        fitFunc = @(x) detrendOut.a*exp(detrendOut.b*x);

    % Fill in the time vector to get only the transient behaviour of the system
        newy = fitFunc(timeData)'+DCgain;

    % Plot the fit and the compensated data for the user to check
        % figure(baseFig_sim+999);clf;hold on
        % set(gcf,'position',[100,100,1000,350])
        %     subplot(121);hold on
        %         plot(timeData,newy);grid minor
        %         plot(timeData,trainDataFilt.y,'--')
        %             xlabel('Time [s]')
        %             ylabel('Temperature [degC]')
        %             title(['Exponential fit on simulated temperature data'])
        %             legend('Exponential fit','Simulation data','location','best')
        %     subplot(122);hold on
        %         plot(timeData,newy-trainDataFilt.y);grid minor
        %             xlabel('Time [s]')
        %             ylabel('Temperature fit error [degC]')
        %             title(['Exponential fit error on simulated temperature data'])
        %             legend(['RMSE: ',num2str(rms(newy-trainDataFilt.y)),'degC'])
        % print('./_Figures/OL/exponentialFit','-depsc')

    % Plot the fit and the compensated data for the user to check (seperate plots)
    %     figure(3111);clf;hold on
    %     set(gcf,'position',[100,100,500,200])
    %         plot(timeData,newy);grid minor
    %         plot(timeData,trainDataFilt.y,'--')
    %             xlabel('Time [s]')
    %             ylabel('Temperature [degC]')
    %             title(['Exponential fit on simulated temperature data'])
    %             legend('Exponential fit','Simulation data','location','best')
    %     print('./_Figures/OL/exponentialFit_OL','-depsc')
    %     figure(3112);clf;hold on
    %     set(gcf,'position',[100,100,500,200])
    %         plot(timeData,newy-trainDataFilt.y);grid minor
    %             xlabel('Time [s]')
    %             ylabel('Temperature fit error [degC]')
    %             title(['Exponential fit error on simulated temperature data'])
    %             legend(['RMSE: ',num2str(rms(newy-trainDataFilt.y)),'degC'])
    %     print('./_Figures/OL/exponentialFit_diff','-depsc')

    % Use the data for new datasets
        sysID.OLdata.train.id = iddata(trainData.y-newy+DCgain,trainData.u,Ts);
        sysID.OLdata.train.y  = trainData.y-newy+DCgain;
        sysID.OLdata.train.u  = trainData.u;
        
        sysID.OLdataFilt.train.id = iddata(trainDataFilt.y-newy+DCgain,trainDataFilt.u,Ts);
        sysID.OLdataFilt.train.y  = trainDataFilt.y-newy+DCgain;
        sysID.OLdataFilt.train.u  = trainDataFilt.u;

    % Using the fit to make a transfer function
        syms t
        F = laplace(detrendOut.a/mean(trainDataFilt.u)*(1-exp(detrendOut.b*t)));
        [num,den] = numden(F);
        num = sym2poly(num);
        den = sym2poly(den);
        sysID.par.OL.time.firstApproxExp = tf(num,den(1:2));

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% From this point on, the data with the transient removed is used %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% FRF measurements using the new state vector and transient compensated data
% Make the FRF's using the raw compensated data
    [sysID.nonPar.OL.frfComp.trd.raw, ~]  = makeOpenLoopFRF_sysIdent(sysID.OLdata.train.y, sysID.OLdata.train.u, fs);

% Make the FRF's using the filtered compensated data
    [sysID.nonPar.OL.frfComp.trd.filt, ~] = makeOpenLoopFRF_sysIdent(sysID.OLdataFilt.train.y, sysID.OLdata.train.u, fs);

% Define and do the non-parametric identification using the LPM
    polyOrder = 6; % The order of the polynomial that is fitted to the data
    locality  = 8; % Amount of points (pos & neg) to consider around apprx freq
    
    disp('Starting OL LPM identification')
        [sysID.nonPar.OL.frfComp.lpm.raw,~]  = sysID_LPM(sysID.OLdata.train.y     ,sysID.OLdata.train.u'     ,fs,locality,polyOrder);
        [sysID.nonPar.OL.frfComp.lpm.filt,~] = sysID_LPM(sysID.OLdataFilt.train.y,sysID.OLdata.train.u',fs,locality,polyOrder);
    disp('Finished OL LPM identification')

%% Parametric system identification, using time data, pre-requisites

% Definitions for identification
    nx = 3; % Model order for fixed order identification

% nx order initial system
    initA = sysID.par.OL.time.firstApprox.raw.A.*eye(nx);
    initB = zeros(nx,1);initB(1) = sysID.par.OL.time.firstApprox.raw.B;
    initC = zeros(1,nx);initC(1) = sysID.par.OL.time.firstApprox.raw.C;
    
    init_sys = idss(initA,initB,[1 0 0],0,zeros(nx,1),zeros(nx,1),0);
    init_sys.Structure.B.Free = ones(nx,1);
    init_sys.Structure.C.Free = ones(1,nx);
    init_sys.Structure.K.Free = zeros(nx,1);

%% Parametric system identification, using time data

optSS_sim = ssestOptions('Focus','Simulation','EnforceStability',true);
optSS_prd = ssestOptions('Focus','Prediction','EnforceStability',true);

% % All tested parametric identification options using OL data, using simulation focus
%         disp('Open Loop identification using an initial system: in progress')
%     sysID.par.OL.time.initSys.sim.raw  = ssest(sysID.OLdata.train.id,init_sys,optSS_sim);
%     sysID.par.OL.time.initSys.sim.filt = ssest(sysID.OLdataFilt.train.id,init_sys,optSS_sim);
%         disp('Open Loop identification using an initial system: done')
% 
%         disp(['Open Loop identification using a fixed order (nx=',num2str(nx),'): in progress'])
%     sysID.par.OL.time.fixedOrder1.sim.raw  = ssest(sysID.OLdata.train.id,1,optSS_sim);
%     sysID.par.OL.time.fixedOrder1.sim.filt = ssest(sysID.OLdataFilt.train.id,1,optSS_sim);
%     sysID.par.OL.time.fixedOrder2.sim.raw  = ssest(sysID.OLdata.train.id,2,optSS_sim);
%     sysID.par.OL.time.fixedOrder2.sim.filt = ssest(sysID.OLdataFilt.train.id,2,optSS_sim);
%     sysID.par.OL.time.fixedOrder.sim.raw  = ssest(sysID.OLdata.train.id,3,optSS_sim);
%     sysID.par.OL.time.fixedOrder.sim.filt = ssest(sysID.OLdataFilt.train.id,3,optSS_sim);
%     sysID.par.OL.time.fixedOrder4.sim.raw  = ssest(sysID.OLdata.train.id,4,optSS_sim);
%     sysID.par.OL.time.fixedOrder4.sim.filt = ssest(sysID.OLdataFilt.train.id,4,optSS_sim);
%     sysID.par.OL.time.fixedOrder5.sim.raw  = ssest(sysID.OLdata.train.id,5,optSS_sim);
%     sysID.par.OL.time.fixedOrder5.sim.filt = ssest(sysID.OLdataFilt.train.id,5,optSS_sim);
%         disp(['Open Loop identification using a fixed order (nx=',num2str(nx),'): done'])

% All tested parametric identification options using OL data, using prediction focus
        disp('Open Loop identification using an initial system: in progress')
    sysID.par.OL.time.initSys.prd.raw  = ssest(sysID.OLdata.train.id,init_sys,optSS_prd);
    sysID.par.OL.time.initSys.prd.filt = ssest(sysID.OLdataFilt.train.id,init_sys,optSS_prd);
        disp('Open Loop identification using an initial system: done')

        disp(['Open Loop identification using a fixed order (nx=',num2str(nx),'): in progress'])
    sysID.par.OL.time.fixedOrder1.prd.raw  = ssest(sysID.OLdata.train.id,1,optSS_prd);
    sysID.par.OL.time.fixedOrder1.prd.filt = ssest(sysID.OLdataFilt.train.id,1,optSS_prd);
    sysID.par.OL.time.fixedOrder2.prd.raw  = ssest(sysID.OLdata.train.id,2,optSS_prd);
    sysID.par.OL.time.fixedOrder2.prd.filt = ssest(sysID.OLdataFilt.train.id,2,optSS_prd);
    sysID.par.OL.time.fixedOrder.prd.raw  = ssest(sysID.OLdata.train.id,3,optSS_prd);
    sysID.par.OL.time.fixedOrder.prd.filt = ssest(sysID.OLdataFilt.train.id,3,optSS_prd);
    sysID.par.OL.time.fixedOrder4.prd.raw  = ssest(sysID.OLdata.train.id,4,optSS_prd);
    sysID.par.OL.time.fixedOrder4.prd.filt = ssest(sysID.OLdataFilt.train.id,4,optSS_prd);
    sysID.par.OL.time.fixedOrder5.prd.raw  = ssest(sysID.OLdata.train.id,5,optSS_prd);
    sysID.par.OL.time.fixedOrder5.prd.filt = ssest(sysID.OLdataFilt.train.id,5,optSS_prd);
        disp(['Open Loop identification using a fixed order (nx=',num2str(nx),'): done'])

% Parametric identification using straight OL data, using prediction focus only (default)
        disp('Open Loop identification using the raw data: in progress')
    sysID.par.OL.time.straight.prd.raw  = ssest(sysID.OLdata.train.id);
    sysID.par.OL.time.straight.prd.filt = ssest(sysID.OLdataFilt.train.id);
        disp('Open Loop identification using the raw data: done')

%% Validate and compare the identifications in time domain, using x-step ahead prediction and a residual test
%% I think this entire section can go away
% 
% % Choose which dataset to use for the validation
%     % dataSet = sysID.OLdata.train.id;
%     dataSet = sysID.OLdata.valid.id;
%     % dataSet = sysID.OLdata.trans.id;
% 
%     % dataSet_filt = sysID.OLdataFilt.train.id;
%     dataSet_filt = sysID.OLdataFilt.valid.id;
%     % dataSet_filt = sysID.OLdataFilt.trans.id;
% 
% % Choose the amount of steps to use for the prediction (inf=simulation)
%     xStep = 1;
%     % xStep = inf;
% 
% if plotAll
% % Run the predictions and plot them
%     % test = initializedVal_OL(sysID.par.OL.time.initSys.prd.raw,sysID.OLdata.train.id);
%     OLsim_sysID_predTests_V3
% 
% % Run the residual tests and plot them
%     % initializedVal_OL(sysID.par.OL.time.initSys.prd.raw,sysID.OLdata.train.id)
%     OLsim_sysID_resTests_V3
% end

%% Frequency domain identification

% Define how may of the lowest frequency points to remove from the responsedata
Nrem = 0;

% Make the FRF measurements usefull for parametric identification
%   Only the compensated FRFs are used
OL_freqDat_trd      = idfrd(squeeze(sysID.nonPar.OL.frfComp.trd.raw.ResponseData(:,:,Nrem+1:end)) ,sysID.nonPar.OL.frfComp.trd.raw.Frequency(Nrem+1:end,:) ,1);
OL_freqDat_lpm      = idfrd(squeeze(sysID.nonPar.OL.frfComp.lpm.raw.ResponseData(:,:,Nrem+1:end)) ,sysID.nonPar.OL.frfComp.lpm.raw.Frequency(Nrem+1:end,:) ,1);
OL_freqDat_trd_filt = idfrd(squeeze(sysID.nonPar.OL.frfComp.trd.filt.ResponseData(:,:,Nrem+1:end)),sysID.nonPar.OL.frfComp.trd.filt.Frequency(Nrem+1:end,:),1);
OL_freqDat_lpm_filt = idfrd(squeeze(sysID.nonPar.OL.frfComp.lpm.filt.ResponseData(:,:,Nrem+1:end)),sysID.nonPar.OL.frfComp.lpm.filt.Frequency(Nrem+1:end,:),1);

% Perform the identifications
sysID.par.OL.freq.trd.prd.raw  = ssest(OL_freqDat_trd);
sysID.par.OL.freq.lpm.prd.raw  = ssest(OL_freqDat_lpm);

sysID.par.OL.freq.trd.prd.filt = ssest(OL_freqDat_trd_filt);
sysID.par.OL.freq.lpm.prd.filt = ssest(OL_freqDat_lpm_filt);

%% Optional identifications
%   These identifications can show the difference that the model order
%   makes to the fit quality, both on raw and filtered data
% sysID.par.OL.freq.trd1.prd.raw  = ssest(OL_freqDat_trd      ,1 ,optSS_prd);
% sysID.par.OL.freq.lpm1.prd.raw  = ssest(OL_freqDat_lpm      ,1 ,optSS_prd);
% 
% sysID.par.OL.freq.trd2.prd.raw  = ssest(OL_freqDat_trd      ,2 ,optSS_prd);
% sysID.par.OL.freq.lpm2.prd.raw  = ssest(OL_freqDat_lpm      ,2 ,optSS_prd);
% 
% sysID.par.OL.freq.trd3.prd.raw  = ssest(OL_freqDat_trd      ,3 ,optSS_prd);
% sysID.par.OL.freq.lpm3.prd.raw  = ssest(OL_freqDat_lpm      ,3 ,optSS_prd);
% 
% sysID.par.OL.freq.trd4.prd.raw  = ssest(OL_freqDat_trd      ,4 ,optSS_prd);
% sysID.par.OL.freq.lpm4.prd.raw  = ssest(OL_freqDat_lpm      ,4 ,optSS_prd);
% 
% sysID.par.OL.freq.trd5.prd.raw  = ssest(OL_freqDat_trd      ,5 ,optSS_prd);
% sysID.par.OL.freq.lpm5.prd.raw  = ssest(OL_freqDat_lpm      ,5 ,optSS_prd);
% 
% 
% sysID.par.OL.freq.trd1.prd.filt  = ssest(OL_freqDat_trd_filt ,1 ,optSS_prd);
% sysID.par.OL.freq.lpm1.prd.filt  = ssest(OL_freqDat_lpm_filt ,1 ,optSS_prd);
% 
% sysID.par.OL.freq.trd2.prd.filt  = ssest(OL_freqDat_trd_filt ,2 ,optSS_prd);
% sysID.par.OL.freq.lpm2.prd.filt  = ssest(OL_freqDat_lpm_filt ,2 ,optSS_prd);
% 
% sysID.par.OL.freq.trd3.prd.filt  = ssest(OL_freqDat_trd_filt ,3 ,optSS_prd);
% sysID.par.OL.freq.lpm3.prd.filt  = ssest(OL_freqDat_lpm_filt ,3 ,optSS_prd);
% 
% sysID.par.OL.freq.trd4.prd.filt  = ssest(OL_freqDat_trd_filt ,4 ,optSS_prd);
% sysID.par.OL.freq.lpm4.prd.filt  = ssest(OL_freqDat_lpm_filt ,4 ,optSS_prd);
% 
% sysID.par.OL.freq.trd5.prd.filt  = ssest(OL_freqDat_trd_filt ,5 ,optSS_prd);
% sysID.par.OL.freq.lpm5.prd.filt  = ssest(OL_freqDat_lpm_filt ,5 ,optSS_prd);


