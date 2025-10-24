
Range = [1e-3 0.1];

figure(baseFig_sim+201);clf
    set(gcf,'position',[500 100 900 500])
    sgtitle(['Comparison between original (filtered) data and transient compensated (filtered) data, ',num2str(tMeas),' hours'])
        subplot(221)
            simpleBodemag(sysID.nonPar.OL.frf.lpm.raw,'Hz',lw);hold on;grid minor
            simpleBodemag(sysID.nonPar.OL.frfComp.lpm.raw,'Hz',lw)
            simpleBodemag(G(1,1),'Hz',lw,'g--')
            simpleBodemag(G(1,1)-sysID.nonPar.OL.frf.lpm.raw,'Hz',0.5,'b--')
            simpleBodemag(G(1,1)-sysID.nonPar.OL.frfComp.lpm.raw,'Hz',0.5,'r--')
            xlim(Range)
                title(['LPM based FRF'])
                legend('Original data','Transient compensated data','Model','location','best')
        subplot(223)
            simpleBodephase(sysID.nonPar.OL.frf.lpm.raw,'Hz',lw,'wrap');hold on;grid minor
            simpleBodephase(sysID.nonPar.OL.frfComp.lpm.raw,'Hz',lw,'wrap')
            simpleBodephase(G(1,1)-sysID.nonPar.OL.frf.lpm.raw,'Hz',0.5,'wrap','b--')
            simpleBodephase(G(1,1)-sysID.nonPar.OL.frfComp.lpm.raw,'Hz',0.5,'wrap','r--')
            simpleBodephase(G(1,1),'Hz',lw,'g--','wrap');
            xlim(Range)
                legend('Original data','Transient compensated data','Model','location','northeast')
        subplot(222)
            simpleBodemag(sysID.nonPar.OL.frf.trd.raw,'Hz',lw);hold on;grid minor
            simpleBodemag(sysID.nonPar.OL.frfComp.trd.raw,'Hz',lw)
            simpleBodemag(G(1,1)-sysID.nonPar.OL.frf.trd.raw,'Hz',0.5,'b--')
            simpleBodemag(G(1,1)-sysID.nonPar.OL.frfComp.trd.raw,'Hz',0.5,'r--')
            simpleBodemag(G(1,1),'Hz',lw,'g--')
            xlim(Range)
                title(['Traditional method based FRF'])
                legend('Original data','Transient compensated data','Model','location','best')
        subplot(224)
            simpleBodephase(sysID.nonPar.OL.frf.trd.raw,'Hz',lw,'wrap');hold on;grid minor
            simpleBodephase(sysID.nonPar.OL.frfComp.trd.raw,'Hz',lw,'wrap')
            simpleBodephase(G(1,1)-sysID.nonPar.OL.frf.trd.raw,'Hz',0.5,'wrap','b--')
            simpleBodephase(G(1,1)-sysID.nonPar.OL.frfComp.trd.raw,'Hz',0.5,'wrap','r--')
            simpleBodephase(G(1,1),'Hz',lw,'g--','wrap')
            xlim(Range)
                legend('Original data','Transient compensated data','Model','location','northeast')
    print('./_Figures/OL/nonPar/orgCompensated_FRF_1','-depsc')
%%
figure(baseFig_sim+202);clf
    set(gcf,'position',[500 100 900 500])
    sgtitle(['Comparison between original (filtered) data and transient compensated (filtered) data, ',num2str(tMeas),' hours'])
        subplot(221)
            xlim([1e-4 0.1])
                title(['LPM based FRF'])
                legend('Original data','Transient compensated data','Model','location','southwest')
        subplot(223)
            simpleBodephase(G(1,1)-sysID.nonPar.OL.frf.lpm.raw,'Hz',lw,'wrap');hold on;grid minor
            simpleBodephase(G(1,1)-sysID.nonPar.OL.frfComp.lpm.raw,'Hz',lw,'wrap')
            xlim([1e-4 0.1])
                legend('Original data','Transient compensated data','Model','location','southwest')
        subplot(222)
            simpleBodemag(G(1,1)-sysID.nonPar.OL.frf.trd.raw,'Hz',lw);hold on;grid minor
            simpleBodemag(G(1,1)-sysID.nonPar.OL.frfComp.trd.raw,'Hz',lw)
            xlim([1e-4 0.1])
                title(['Traditional method based FRF'])
                legend('Original data','Transient compensated data','Model','location','southwest')
        subplot(224)
            simpleBodephase(G(1,1)-sysID.nonPar.OL.frf.trd.raw,'Hz',lw,'wrap');hold on;grid minor
            simpleBodephase(G(1,1)-sysID.nonPar.OL.frfComp.trd.raw,'Hz',lw,'wrap')
            xlim([1e-4 0.1])
                legend('Original data','Transient compensated data','Model','location','southwest')
    print('./_Figures/OL/nonPar/orgCompensated_FRF_1','-depsc')









