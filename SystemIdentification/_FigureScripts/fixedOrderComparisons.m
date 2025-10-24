
%%
figure(baseFig_sim+2);clf
set(gcf,'Position',[100 100 500 450])
    subplot(211);
        simpleBodemag(sysID.par.OL.time.fixedOrder1.prd.filt,'Hz',lw,bodeRange);hold on;grid minor
        simpleBodemag(sysID.par.OL.time.fixedOrder2.prd.filt,'Hz',lw,bodeRange)
        simpleBodemag(sysID.par.OL.time.fixedOrder.prd.filt ,'Hz',lw,bodeRange)
        simpleBodemag(sysID.par.OL.time.fixedOrder4.prd.filt,'Hz',lw,bodeRange)
        simpleBodemag(sysID.par.OL.time.fixedOrder5.prd.filt,'Hz',lw,bodeRange)
        simpleBodemag(G(1,1),'Hz','g--',lw,bodeRange)
            title('Results from different fixed order identifications')
            legend('nx=1','nx=2','nx=3','nx=4','nx=5','location','southwest')
    subplot(212);
        simpleBodemag(G(1,1)-sysID.par.OL.time.fixedOrder1.prd.filt,'Hz',lw,bodeRange);hold on;grid minor
        simpleBodemag(G(1,1)-sysID.par.OL.time.fixedOrder2.prd.filt,'Hz',lw,bodeRange)
        simpleBodemag(G(1,1)-sysID.par.OL.time.fixedOrder.prd.filt ,'Hz',lw,bodeRange)
        simpleBodemag(G(1,1)-sysID.par.OL.time.fixedOrder4.prd.filt,'Hz',lw,bodeRange)
        simpleBodemag(G(1,1)-sysID.par.OL.time.fixedOrder5.prd.filt,'Hz',lw,bodeRange)
            ylabel('Magnitude error [dB]')
            title('Magnitude error of different fixed order identifications')
            legend('nx=1','nx=2','nx=3','nx=4','nx=5','location','southwest')
print('./_Figures/OL/parTime/fixedOrderComparison_bodes','-depsc')
    
preds_fixedOrders = compare(sysID.OLdataFilt.valid.id,sysID.par.OL.time.fixedOrder1.sim.filt,sysID.par.OL.time.fixedOrder2.sim.filt,sysID.par.OL.time.fixedOrder.sim.filt,sysID.par.OL.time.fixedOrder4.sim.filt,sysID.par.OL.time.fixedOrder5.sim.filt,1);
sims_fixedOrders  = compare(sysID.OLdataFilt.valid.id,sysID.par.OL.time.fixedOrder1.sim.filt,sysID.par.OL.time.fixedOrder2.sim.filt,sysID.par.OL.time.fixedOrder.sim.filt,sysID.par.OL.time.fixedOrder4.sim.filt,sysID.par.OL.time.fixedOrder5.sim.filt);

preds_fixedOrders_err{1} = sysID.OLdataFilt.valid.y-squeeze(preds_fixedOrders{1}.y);
preds_fixedOrders_err{2} = sysID.OLdataFilt.valid.y-squeeze(preds_fixedOrders{2}.y);
preds_fixedOrders_err{3} = sysID.OLdataFilt.valid.y-squeeze(preds_fixedOrders{3}.y);
preds_fixedOrders_err{4} = sysID.OLdataFilt.valid.y-squeeze(preds_fixedOrders{4}.y);
preds_fixedOrders_err{5} = sysID.OLdataFilt.valid.y-squeeze(preds_fixedOrders{5}.y);

sims_fixedOrders_err{1} = sysID.OLdataFilt.valid.y-squeeze(sims_fixedOrders{1}.y);
sims_fixedOrders_err{2} = sysID.OLdataFilt.valid.y-squeeze(sims_fixedOrders{2}.y);
sims_fixedOrders_err{3} = sysID.OLdataFilt.valid.y-squeeze(sims_fixedOrders{3}.y);
sims_fixedOrders_err{4} = sysID.OLdataFilt.valid.y-squeeze(sims_fixedOrders{4}.y);
sims_fixedOrders_err{5} = sysID.OLdataFilt.valid.y-squeeze(sims_fixedOrders{5}.y);

figure(baseFig_sim+3);clf
set(gcf,'Position',[100 100 1000 450])
    subplot(221);hold on;grid minor
        plot(sysID.OLdata.valid.tVec,sysID.OLdataFilt.valid.y)
        plot(sysID.OLdata.valid.tVec,squeeze(preds_fixedOrders{1}.y))
        plot(sysID.OLdata.valid.tVec,squeeze(preds_fixedOrders{2}.y))
        plot(sysID.OLdata.valid.tVec,squeeze(preds_fixedOrders{3}.y))
        plot(sysID.OLdata.valid.tVec,squeeze(preds_fixedOrders{4}.y))
        plot(sysID.OLdata.valid.tVec,squeeze(preds_fixedOrders{5}.y))
            xlabel('Time [s]')
            ylabel('Temperature [degC]')
            title('1-Step ahead prediction using several fixed order models')
            legend('Validation data','nx=1','nx=2','nx=3','nx=4','nx=5','location','best')
    subplot(223);hold on;grid minor
        plot(sysID.OLdata.valid.tVec,preds_fixedOrders_err{1})
        plot(sysID.OLdata.valid.tVec,preds_fixedOrders_err{2})
        plot(sysID.OLdata.valid.tVec,preds_fixedOrders_err{3})
        plot(sysID.OLdata.valid.tVec,preds_fixedOrders_err{4})
        plot(sysID.OLdata.valid.tVec,preds_fixedOrders_err{5})
            xlabel('Time [s]')
            ylabel('Temperature error [degC]')
            title('1-Step ahead prediction error using several fixed order models')
            legend(['nx=1, RMSE: ',num2str(rms(preds_fixedOrders_err{1}),3),'degC'],...
                   ['nx=2, RMSE: ',num2str(rms(preds_fixedOrders_err{2}),3),'degC'],...
                   ['nx=3, RMSE: ',num2str(rms(preds_fixedOrders_err{3}),3),'degC'],...
                   ['nx=4, RMSE: ',num2str(rms(preds_fixedOrders_err{4}),3),'degC'],...
                   ['nx=5, RMSE: ',num2str(rms(preds_fixedOrders_err{5}),3),'degC'],...
                   'location','best')
    subplot(222);hold on;grid minor
        plot(sysID.OLdata.valid.tVec,sysID.OLdataFilt.valid.y)
        plot(sysID.OLdata.valid.tVec,squeeze(sims_fixedOrders{1}.y))
        plot(sysID.OLdata.valid.tVec,squeeze(sims_fixedOrders{2}.y))
        plot(sysID.OLdata.valid.tVec,squeeze(sims_fixedOrders{3}.y))
        plot(sysID.OLdata.valid.tVec,squeeze(sims_fixedOrders{4}.y))
        plot(sysID.OLdata.valid.tVec,squeeze(sims_fixedOrders{5}.y))
            xlabel('Time [s]')
            ylabel('Temperature [degC]')
            title('Simulation using several fixed order models')
            legend('Validation data','nx=1','nx=2','nx=3','nx=4','nx=5','location','best')
    subplot(224);hold on;grid minor
        plot(sysID.OLdata.valid.tVec,sims_fixedOrders_err{1})
        plot(sysID.OLdata.valid.tVec,sims_fixedOrders_err{2})
        plot(sysID.OLdata.valid.tVec,sims_fixedOrders_err{3})
        plot(sysID.OLdata.valid.tVec,sims_fixedOrders_err{4})
        plot(sysID.OLdata.valid.tVec,sims_fixedOrders_err{5})
            xlabel('Time [s]')
            ylabel('Temperature error [degC]')
            title('Simulation error using several fixed order models')
            legend(['nx=1, RMSE: ',num2str(rms(sims_fixedOrders_err{1}),3),'degC'],...
                   ['nx=2, RMSE: ',num2str(rms(sims_fixedOrders_err{2}),3),'degC'],...
                   ['nx=3, RMSE: ',num2str(rms(sims_fixedOrders_err{3}),3),'degC'],...
                   ['nx=4, RMSE: ',num2str(rms(sims_fixedOrders_err{4}),3),'degC'],...
                   ['nx=5, RMSE: ',num2str(rms(sims_fixedOrders_err{5}),3),'degC'],...
                   'location','best')
print('./_Figures/OL/parTime/fixedOrderComparison_predSim','-depsc')

