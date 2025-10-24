clc

%% Visualize and compare the identifications in time domain
% dataSet = OL_dat;
% dataSet = OL_dat_val;
% dataSet = OL_dat_trns;

% dataSet_filt = OL_dat_filt;
% dataSet_filt = OL_dat_val_filt;
% dataSet_filt = OL_dat_trns_filt;

% figure(baseFig_sim+401);clf
%     compare(dataSet_filt,sysID.par.initSys.sim.OL.raw,sysID.par.fixedOrder.sim.OL.raw,sysID.par.firstAprrox.OL.raw,sysID.par.straight.sim.OL.raw,xStep);grid minor
%         legend('Real data','Initial system','Fixed order','First order approx','Raw data')
%     print('./_Figures/OL/parTime/compareSimAll_raw_1','-depsc')
% 
% figure(baseFig_sim+402);clf
%     compare(dataSet_filt,sysID.par.initSys.sim.OL.filt,sysID.par.fixedOrder.sim.OL.filt,sysID.par.firstAprrox.OL.filt,sysID.par.straight.sim.OL.filt,xStep);grid minor
%         legend('Real data','Initial system','Fixed order','First order approx','Raw data')
%     print('./_Figures/OL/parTime/compareSimAll_filt_1','-depsc')

figure(baseFig_sim+411);clf
    dataNames = {'Initial system','Fixed order','First order approx','Raw data'};
    compareRMS(dataSet,dataNames,sysID.par.initSys.sim.OL.raw,sysID.par.fixedOrder.sim.OL.raw,sysID.par.firstAprrox.OL.raw,sysID.par.straight.sim.OL.raw,xStep);
    print('./_Figures/OL/parTime/compareSimAll_RMSE_raw_1','-depsc')

figure(baseFig_sim+412);clf
    dataNames = {'Initial system','Fixed order','First order approx','Raw data'};
    compareRMS(dataSet,dataNames,sysID.par.initSys.sim.OL.filt,sysID.par.fixedOrder.sim.OL.filt,sysID.par.firstAprrox.OL.filt,sysID.par.straight.sim.OL.filt,xStep);
    print('./_Figures/OL/parTime/compareSimAll_RMSE_filt_1','-depsc')

%% x-step ahead prediction tests
% xStep = 1;
% xStep = inf;

figure(baseFig_sim+501);clf
    subplot(241)
        compare(dataSet,sysID.par.initSys.sim.OL.raw,sysID.par.initSys.sim.OL.filt,xStep);grid minor
            legend('Real (raw) data','Raw based model','Filtered based model',Location='best')
            title([num2str(xStep),'-step prediction on raw data, initSys'])
    subplot(245)
        compare(dataSet_filt,sysID.par.initSys.sim.OL.raw,sysID.par.initSys.sim.OL.filt,xStep);grid minor
            legend('Real (filtered) data','Raw based model','Filtered based model',Location='best')
            title([num2str(xStep),'-step prediction on filt data, initSys'])
    subplot(242)
        compare(dataSet,sysID.par.fixedOrder.sim.OL.raw,sysID.par.fixedOrder.sim.OL.filt,xStep);grid minor
            legend('Real (raw) data','Raw based model','Filtered based model',Location='best')
            title([num2str(xStep),'-step prediction on raw data, fixed order'])
    subplot(246)
        compare(dataSet_filt,sysID.par.fixedOrder.sim.OL.raw,sysID.par.fixedOrder.sim.OL.filt,xStep);grid minor
            legend('Real (filtered) data','Raw based model','Filtered based model',Location='best')
            title([num2str(xStep),'-step prediction on filt data, fixed order'])
    subplot(243)
        compare(dataSet,sysID.par.firstAprrox.OL.raw,sysID.par.firstAprrox.OL.filt,xStep);grid minor
            legend('Real (raw) data','Raw based model','Filtered based model',Location='best')
            title([num2str(xStep),'-step prediction on raw data, first order approx'])
    subplot(247)
        compare(dataSet_filt,sysID.par.firstAprrox.OL.raw,sysID.par.firstAprrox.OL.filt,xStep);grid minor
            legend('Real (filtered) data','Raw based model','Filtered based model',Location='best')
            title([num2str(xStep),'-step prediction on filt data, first order approx'])
    subplot(244)
        compare(dataSet,sysID.par.straight.sim.OL.raw,sysID.par.straight.sim.OL.filt,xStep);grid minor
            legend('Real (raw) data','Raw based model','Filtered based model',Location='best')
            title([num2str(xStep),'-step prediction on raw data, straight data'])
    subplot(248)
        compare(dataSet_filt,sysID.par.straight.sim.OL.raw,sysID.par.straight.sim.OL.filt,xStep);grid minor
            legend('Real (filtered) data','Raw based model','Filtered based model',Location='best')
            title([num2str(xStep),'-step prediction on filt data, straight data'])
    print('./_Figures/OL/parTime/compareSimAll_1','-depsc')

figure(baseFig_sim+511);clf
    subplot(241)
        dataNames = {'Raw based model','Filtered based model'};
        compareRMS(dataSet,dataNames,sysID.par.initSys.sim.OL.raw,sysID.par.initSys.sim.OL.filt,xStep);
            title([num2str(xStep),'-step prediction on raw data, initSys'])
    subplot(245)
        dataNames = {'Raw based model','Filtered based model'};
        compareRMS(dataSet_filt,dataNames,sysID.par.initSys.sim.OL.raw,sysID.par.initSys.sim.OL.filt,xStep);
            title([num2str(xStep),'-step prediction on filt data, initSys'])
    subplot(242)
        dataNames = {'Raw based model','Filtered based model'};
        compareRMS(dataSet,dataNames,sysID.par.fixedOrder.sim.OL.raw,sysID.par.fixedOrder.sim.OL.filt,xStep);
            title([num2str(xStep),'-step prediction on raw data, fixed order'])
    subplot(246)
        dataNames = {'Raw based model','Filtered based model'};
        compareRMS(dataSet_filt,dataNames,sysID.par.fixedOrder.sim.OL.raw,sysID.par.fixedOrder.sim.OL.filt,xStep);
            title([num2str(xStep),'-step prediction on filt data, fixed order'])
    subplot(243)
        dataNames = {'Raw based model','Filtered based model'};
        compareRMS(dataSet,dataNames,sysID.par.firstAprrox.OL.raw,sysID.par.firstAprrox.OL.filt,xStep);
            title([num2str(xStep),'-step prediction on raw data, first order approx'])
    subplot(247)
        dataNames = {'Raw based model','Filtered based model'};
        compareRMS(dataSet_filt,dataNames,sysID.par.firstAprrox.OL.raw,sysID.par.firstAprrox.OL.filt,xStep);
            title([num2str(xStep),'-step prediction on filt data, first order approx'])
    subplot(244)
        dataNames = {'Raw based model','Filtered based model'};
        compareRMS(dataSet,dataNames,sysID.par.straight.sim.OL.raw,sysID.par.straight.sim.OL.filt,xStep);
            title([num2str(xStep),'-step prediction on raw data, straight data'])
    subplot(248)
        dataNames = {'Raw based model','Filtered based model'};
        compareRMS(dataSet_filt,dataNames,sysID.par.straight.sim.OL.raw,sysID.par.straight.sim.OL.filt,xStep);
            title([num2str(xStep),'-step prediction on filt data, straight data'])
    print('./_Figures/OL/parTime/compareSimAll_RMSE_1','-depsc')









