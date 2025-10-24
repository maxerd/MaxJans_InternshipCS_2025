
%%

figure(baseFig_sim+4);clf
set(gcf,'Position',[100 100 500 450])
    subplot(211);
        simpleBodemag(sysID.par.OL.time.firstApprox.filt    ,'Hz',lw,bodeRange);hold on;grid minor
        simpleBodemag(sysID.par.OL.time.fixedOrder1.prd.filt,'Hz',lw,bodeRange)
        simpleBodemag(sysID.par.OL.time.fixedOrder5.prd.filt,'Hz',lw,bodeRange)
        simpleBodemag(sysID.par.OL.time.initSys.prd.filt    ,'Hz',lw,bodeRange)
        simpleBodemag(sysID.par.OL.time.straight.prd.filt   ,'Hz',lw,bodeRange)
        simpleBodemag(G(1,1),'Hz','g--',lw,bodeRange)
            title('Results from different identification approaches')
            legend('First order approx.','Fixed order, nx=1','Fixed order, nx=5','Initialized system','Straight (filtered) data','Model','location','southwest')
    subplot(212);
        simpleBodemag(G(1,1)-sysID.par.OL.time.firstApprox.filt    ,'Hz',lw,bodeRange);hold on;grid minor
        simpleBodemag(G(1,1)-sysID.par.OL.time.fixedOrder1.prd.filt,'Hz',lw,bodeRange)
        simpleBodemag(G(1,1)-sysID.par.OL.time.fixedOrder5.prd.filt,'Hz',lw,bodeRange)
        simpleBodemag(G(1,1)-sysID.par.OL.time.initSys.prd.filt    ,'Hz',lw,bodeRange)
        simpleBodemag(G(1,1)-sysID.par.OL.time.straight.prd.filt   ,'Hz',lw,bodeRange)
            ylabel('Magnitude error [dB]')
            title('Magnitude error of identification approaches compared to model')
            legend('First order approx.','Fixed order, nx=1','Fixed order, nx=5','Initialized system','Straight (filtered) data','Model','location','south')
print('./_Figures/OL/parTime/identificationsComparison_bodes','-depsc')
    

preds_timeApproaches = compare(sysID.OLdataFilt.valid.id,...
                            sysID.par.OL.time.firstApprox.filt,...
                            sysID.par.OL.time.fixedOrder1.prd.filt,...
                            sysID.par.OL.time.fixedOrder5.prd.filt,...
                            sysID.par.OL.time.initSys.prd.filt,...
                            sysID.par.OL.time.straight.prd.filt,1);
sims_timeApproaches  = compare(sysID.OLdataFilt.valid.id,...
                            sysID.par.OL.time.firstApprox.filt,...
                            sysID.par.OL.time.fixedOrder1.prd.filt,...
                            sysID.par.OL.time.fixedOrder5.prd.filt,...
                            sysID.par.OL.time.initSys.prd.filt,...
                            sysID.par.OL.time.straight.prd.filt);

preds_timeApproaches_err{1} = sysID.OLdataFilt.valid.y-squeeze(preds_timeApproaches{1}.y);
preds_timeApproaches_err{2} = sysID.OLdataFilt.valid.y-squeeze(preds_timeApproaches{2}.y);
preds_timeApproaches_err{3} = sysID.OLdataFilt.valid.y-squeeze(preds_timeApproaches{3}.y);
preds_timeApproaches_err{4} = sysID.OLdataFilt.valid.y-squeeze(preds_timeApproaches{4}.y);
preds_timeApproaches_err{5} = sysID.OLdataFilt.valid.y-squeeze(preds_timeApproaches{5}.y);

sims_timeApproaches_err{1} = sysID.OLdataFilt.valid.y-squeeze(sims_timeApproaches{1}.y);
sims_timeApproaches_err{2} = sysID.OLdataFilt.valid.y-squeeze(sims_timeApproaches{2}.y);
sims_timeApproaches_err{3} = sysID.OLdataFilt.valid.y-squeeze(sims_timeApproaches{3}.y);
sims_timeApproaches_err{4} = sysID.OLdataFilt.valid.y-squeeze(sims_timeApproaches{4}.y);
sims_timeApproaches_err{5} = sysID.OLdataFilt.valid.y-squeeze(sims_timeApproaches{5}.y);

figure(baseFig_sim+5);clf
set(gcf,'Position',[100 100 1200 650])
    subplot(221);hold on;grid minor
        plot(sysID.OLdata.valid.tVec,sysID.OLdataFilt.valid.y)
        plot(sysID.OLdata.valid.tVec,squeeze(preds_timeApproaches{1}.y))
        plot(sysID.OLdata.valid.tVec,squeeze(preds_timeApproaches{2}.y))
        plot(sysID.OLdata.valid.tVec,squeeze(preds_timeApproaches{3}.y))
        plot(sysID.OLdata.valid.tVec,squeeze(preds_timeApproaches{4}.y))
        plot(sysID.OLdata.valid.tVec,squeeze(preds_timeApproaches{5}.y))
            xlabel('Time [s]')
            ylabel('Temperature [degC]')
            title('1-Step ahead prediction using several fixed order models')
            legend('Validation data','First order approx.','Fixed order, nx=1','Fixed order, nx=5','Initialized system','Straight (filtered) data','location','best')
    subplot(223);hold on;grid minor
        plot(sysID.OLdata.valid.tVec,preds_timeApproaches_err{1})
        plot(sysID.OLdata.valid.tVec,preds_timeApproaches_err{2})
        plot(sysID.OLdata.valid.tVec,preds_timeApproaches_err{3})
        plot(sysID.OLdata.valid.tVec,preds_timeApproaches_err{4})
        plot(sysID.OLdata.valid.tVec,preds_timeApproaches_err{5})
            xlabel('Time [s]')
            ylabel('Temperature error [degC]')
            title('1-Step ahead prediction error using several fixed order models')
            legend('First order approx.','Fixed order, nx=1','Fixed order, nx=5','Initialized system','Straight (filtered) data')
    subplot(222);hold on;grid minor
        plot(sysID.OLdata.valid.tVec,sysID.OLdataFilt.valid.y)
        plot(sysID.OLdata.valid.tVec,squeeze(sims_timeApproaches{1}.y))
        plot(sysID.OLdata.valid.tVec,squeeze(sims_timeApproaches{2}.y))
        plot(sysID.OLdata.valid.tVec,squeeze(sims_timeApproaches{3}.y))
        plot(sysID.OLdata.valid.tVec,squeeze(sims_timeApproaches{4}.y))
        plot(sysID.OLdata.valid.tVec,squeeze(sims_timeApproaches{5}.y))
            xlabel('Time [s]')
            ylabel('Temperature [degC]')
            title('Simulation using several fixed order models')
            legend('Validation data','First order approx.','Fixed order, nx=1','Fixed order, nx=5','Initialized system','Straight (filtered) data','location','best')
    subplot(224);hold on;grid minor
        plot(sysID.OLdata.valid.tVec,sims_timeApproaches_err{1})
        plot(sysID.OLdata.valid.tVec,sims_timeApproaches_err{2})
        plot(sysID.OLdata.valid.tVec,sims_timeApproaches_err{3})
        plot(sysID.OLdata.valid.tVec,sims_timeApproaches_err{4})
        plot(sysID.OLdata.valid.tVec,sims_timeApproaches_err{5})
            xlabel('Time [s]')
            ylabel('Temperature error [degC]')
            title('Simulation error using several fixed order models')
            legend('First order approx.','Fixed order, nx=1','Fixed order, nx=5','Initialized system','Straight (filtered) data')
print('./_Figures/OL/parTime/identificationsComparison_predSim','-depsc')

disp(' ')
disp('RMS prediction errors of time based identification approaches:')
disp(['First order approx., RMSE: ',     num2str(rms(preds_timeApproaches_err{1}),3),'degC'])
disp(['Fixed order, nx=1, RMSE: ',       num2str(rms(preds_timeApproaches_err{2}),3),'degC'])
disp(['Fixed order, nx=5, RMSE: ',       num2str(rms(preds_timeApproaches_err{3}),3),'degC'])
disp(['Initialized system, RMSE: ',      num2str(rms(preds_timeApproaches_err{4}),3),'degC'])
disp(['Straight (filtered) data, RMSE: ',num2str(rms(preds_timeApproaches_err{5}),3),'degC'])
disp(' ')
disp('RMS simulation errors of time based identification approaches:')
disp(['First order approx., RMSE: ',     num2str(rms(sims_timeApproaches_err{1}),3),'degC'])
disp(['Fixed order, nx=1, RMSE: ',       num2str(rms(sims_timeApproaches_err{2}),3),'degC'])
disp(['Fixed order, nx=5, RMSE: ',       num2str(rms(sims_timeApproaches_err{3}),3),'degC'])
disp(['Initialized system, RMSE: ',      num2str(rms(sims_timeApproaches_err{4}),3),'degC'])
disp(['Straight (filtered) data, RMSE: ',num2str(rms(sims_timeApproaches_err{5}),3),'degC'])

figure(baseFig_sim+6);clf
set(gcf,'Position',[100 100 800 450])
    resid(sysID.OLdataFilt.valid.id,...
          sysID.par.OL.time.firstApprox.filt,...
          sysID.par.OL.time.fixedOrder1.prd.filt,...
          sysID.par.OL.time.fixedOrder5.prd.filt,...
          sysID.par.OL.time.initSys.prd.filt,...
          sysID.par.OL.time.straight.prd.filt);grid minor
        legend('First order approx.','Fixed order, nx=1','Fixed order, nx=5','Initialized system','Straight (filtered) data')
print('./_Figures/OL/parTime/identificationsComparison_resid','-depsc')


return
%%

% simpleResid(sysID.OLdataFilt.valid.id,sysID.par.OL.time.straight.prd.filt,99)
simpleResid(sysID.OLdataFilt.valid.id,sysID.par.OL.time.initSys.prd.filt,99)
% simpleResid(sysID.OLdataFilt.valid.id,sysID.par.OL.time.fixedOrder1.prd.filt,99)
% simpleResid(sysID.OLdataFilt.valid.id,sysID.par.OL.time.firstApprox.filt,99)

figure(3);clf;
    % resid(sysID.OLdataFilt.valid.id,sysID.par.OL.time.straight.prd.filt)
    resid(sysID.OLdataFilt.valid.id,sysID.par.OL.time.initSys.prd.filt)
    % resid(sysID.OLdataFilt.valid.id,sysID.par.OL.time.fixedOrder1.prd.filt)
    % resid(sysID.OLdataFilt.valid.id,sysID.par.OL.time.firstApprox.filt)

function simpleResid(data,model,confBound)
% [vals1,testVal1] = resid(sysID.OLdataFilt.valid.id,sysID.par.OL.time.firstApprox.filt,length(sysID.OLdataFilt.valid.u));
% [vals2,testVal2] = resid(sysID.OLdataFilt.valid.id,sysID.par.OL.time.fixedOrder1.prd.filt,100);
% [vals3,testVal3] = resid(sysID.OLdataFilt.valid.id,sysID.par.OL.time.fixedOrder5.prd.filt,100);
% [vals4,testVal4] = resid(sysID.OLdataFilt.valid.id,sysID.par.OL.time.initSys.prd.filt,100);
[~,testVal5] = resid(data,model,50);


acorr_5_right = ([(testVal5(:,1,1))]);
acorr_5_left = ([testVal5(:,1,1)]);
acorr_5 = ([flip(acorr_5_left(2:end));acorr_5_right]);

xcorr_5_right = ([(testVal5(:,2,1))]);
xcorr_5_left = ([testVal5(:,1,2)]);
xcorr_5 = ([flip(xcorr_5_left(2:end));xcorr_5_right])/sqrt(testVal5(1,1,1)*testVal5(1,2,2));

N  = length(data.y);       % total number of samples
try
    np = length(getpvec(model)); % #params in model
catch
    np = 20;
end
Neff = N - np;                               % effective sample size
if confBound == 99
    conf = 2.576 / sqrt(Neff) / sqrt(testVal5(1,1,1)*testVal5(1,2,2));
    aconf = 2.576 / sqrt(Neff);
elseif confBound == 95
    conf = 1.96 / sqrt(Neff) / sqrt(testVal5(1,1,1)*testVal5(1,2,2));
else
    warning('Confidence bound not known, using 99%')
    conf = 2.576 / sqrt(Neff) / sqrt(testVal5(1,1,1)*testVal5(1,2,2));
end


figure(1);clf;
ax1 = subplot(121);
    stem([-50:50],acorr_5/max(acorr_5));hold on
    xlim([-25 25])
    yline(aconf)
    yline(-aconf)
ax2 = subplot(122);
    stem([-50:50],xcorr_5);hold on
    xlim([-25 25])
    yline(conf)
    yline(-conf)
linkaxes([ax1,ax2],'y')

% figure(2);clf;
%     plot(vals2.y);hold on


end


































