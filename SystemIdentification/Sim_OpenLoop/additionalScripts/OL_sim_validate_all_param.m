
%% Visualize and compare the identifications in frequency domain

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
