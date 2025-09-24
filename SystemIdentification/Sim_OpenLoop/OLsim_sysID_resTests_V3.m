
%% Redisual tests
figure(baseFig_sim+601);clf
    subplot(241)
        resid(dataSet,sysID.par.initSys.sim.OL.raw);grid minor
            % legend('Raw based model','Filtered based model',Location='best')
            title(['Residual tests, initSys, raw'])
    subplot(245)
        resid(dataSet_filt,sysID.par.initSys.sim.OL.raw);grid minor
            % legend('Raw based model','Filtered based model',Location='best')
            title(['Residual tests, initSys, filt'])
    subplot(242)
        resid(dataSet,sysID.par.fixedOrder.sim.OL.raw);grid minor
            % legend('Raw based model','Filtered based model',Location='best')
            title(['Residual tests, fixed order, raw'])
    subplot(246)
        resid(dataSet_filt,sysID.par.fixedOrder.sim.OL.raw);grid minor
            % legend('Raw based model','Filtered based model',Location='best')
            title(['Residual tests, fixed order, filt'])
    subplot(243)
        resid(dataSet,sysID.par.firstAprrox.OL.raw);grid minor
            % legend('Raw based model','Filtered based model',Location='best')
            title(['Residual tests, first order approx, raw'])
    subplot(247)
        resid(dataSet_filt,sysID.par.firstAprrox.OL.raw);grid minor
            % legend('Raw based model','Filtered based model',Location='best')
            title(['Residual tests, first order approx, filt'])
    subplot(244)
        resid(dataSet,sysID.par.straight.sim.OL.raw);grid minor
            % legend('Raw based model','Filtered based model',Location='best')
            title(['Residual tests, straight data, raw'])
    subplot(248)
        resid(dataSet_filt,sysID.par.straight.sim.OL.raw);grid minor
            % legend('Raw based model','Filtered based model',Location='best')
            title(['Residual tests, straight data, filt'])

figure(baseFig_sim+602);clf
    subplot(241)
        resid(dataSet,sysID.par.initSys.sim.OL.filt);grid minor
            % legend('Raw based model','Filtered based model',Location='best')
            title(['Residual tests, initSys, raw'])
    subplot(245)
        resid(dataSet_filt,sysID.par.initSys.sim.OL.filt);grid minor
            % legend('Raw based model','Filtered based model',Location='best')
            title(['Residual tests, initSys, filt'])
    subplot(242)
        resid(dataSet,sysID.par.fixedOrder.sim.OL.filt);grid minor
            % legend('Raw based model','Filtered based model',Location='best')
            title(['Residual tests, fixed order, raw'])
    subplot(246)
        resid(dataSet_filt,sysID.par.fixedOrder.sim.OL.filt);grid minor
            % legend('Raw based model','Filtered based model',Location='best')
            title(['Residual tests, fixed order, filt'])
    subplot(243)
        resid(dataSet,sysID.par.firstAprrox.OL.filt);grid minor
            % legend('Raw based model','Filtered based model',Location='best')
            title(['Residual tests, first order approx, raw'])
    subplot(247)
        resid(dataSet_filt,sysID.par.firstAprrox.OL.filt);grid minor
            % legend('Raw based model','Filtered based model',Location='best')
            title(['Residual tests, first order approx, filt'])
    subplot(244)
        resid(dataSet,sysID.par.straight.sim.OL.filt);grid minor
            % legend('Raw based model','Filtered based model',Location='best')
            title(['Residual tests, straight data, raw'])
    subplot(248)
        resid(dataSet_filt,sysID.par.straight.sim.OL.filt);grid minor
            % legend('Raw based model','Filtered based model',Location='best')
            title(['Residual tests, straight data, filt'])

% figure(baseFig_sim+603);clf
%     subplot(241)
%         resid(dataSet,OLpar_initSys,OLpar_initSys_filt);grid minor
%             legend('Raw based model','Filtered based model',Location='best')
%             title(['Residual tests, initSys, raw'])
%     subplot(245)
%         resid(dataSet_filt,OLpar_initSys,OLpar_initSys_filt);grid minor
%             legend('Raw based model','Filtered based model',Location='best')
%             title(['Residual tests, initSys, filt'])
%     subplot(242)
%         resid(dataSet,OLpar_fixdOrder,OLpar_fixdOrder_filt);grid minor
%             legend('Raw based model','Filtered based model',Location='best')
%             title(['Residual tests, fixed order, raw'])
%     subplot(246)
%         resid(dataSet_filt,OLpar_fixdOrder,OLpar_fixdOrder_filt);grid minor
%             legend('Raw based model','Filtered based model',Location='best')
%             title(['Residual tests, fixed order, filt'])
%     subplot(243)
%         resid(dataSet,OLpar_firstAprrox,OLpar_firstAprrox_filt);grid minor
%             legend('Raw based model','Filtered based model',Location='best')
%             title(['Residual tests, first order approx, raw'])
%     subplot(247)
%         resid(dataSet_filt,OLpar_firstAprrox,OLpar_firstAprrox_filt);grid minor
%             legend('Raw based model','Filtered based model',Location='best')
%             title(['Residual tests, first order approx, filt'])
%     subplot(244)
%         resid(dataSet,OLpar_raw,OLpar_filt);grid minor
%             legend('Raw based model','Filtered based model',Location='best')
%             title(['Residual tests, straight data, raw'])
%     subplot(248)
%         resid(dataSet_filt,OLpar_raw,OLpar_filt);grid minor
%             legend('Raw based model','Filtered based model',Location='best')
%             title(['Residual tests, straight data, filt'])


