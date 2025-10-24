
%% Visualization of the non-parametric identification
% Open loop data
    figure(baseFig_sim+201);clf
    set(gcf,'position',[500 100 900 500])
    sgtitle(['Heater to TM temperature, open loop identification using unfiltered and filtered data, ',num2str(tMeas),' hours'])
        subplot(221)
            simpleBodemag(sysID.nonPar.OL.frf.trd.raw,'Hz',lw);hold on;grid minor
            simpleBodemag(sysID.nonPar.OL.frf.lpm.raw,'Hz',lw);grid minor
            simpleBodemag(G(1,1),'Hz',lw,'g--');grid minor
            xlim([1e-4 0.1])
                title(['Using unfiltered data'])
                legend('Traditional','LPM','Model','location','best')
        subplot(223)
            simpleBodephase(sysID.nonPar.OL.frf.trd.raw,'Hz',lw,'wrap');hold on;grid minor
            simpleBodephase(sysID.nonPar.OL.frf.lpm.raw,'Hz',lw,'wrap');grid minor
            simpleBodephase(G(1,1),'Hz',lw,'g--','wrap');grid minor
            xlim([1e-4 0.1])
                legend('Traditional','LPM','Model','location','best')
        subplot(222)
            simpleBodemag(sysID.nonPar.OL.frf.trd.filt,'Hz',lw);hold on;grid minor
            simpleBodemag(sysID.nonPar.OL.frf.lpm.filt,'Hz',lw);grid minor
            simpleBodemag(G(1,1),'Hz',lw,'g--');grid minor
            xlim([1e-4 0.1])
                title(['Using filtered data'])
                legend('Traditional','LPM','Model','location','best')
        subplot(224)
            simpleBodephase(sysID.nonPar.OL.frf.trd.filt,'Hz',lw,'wrap');hold on;grid minor
            simpleBodephase(sysID.nonPar.OL.frf.lpm.filt,'Hz',lw,'wrap');grid minor
            simpleBodephase(G(1,1),'Hz',lw,'g--','wrap');grid minor
            xlim([1e-4 0.1])
                legend('Traditional','LPM','Model','location','best')
    % print('./_Figures/OL/nonPar/tradAndLPM_FRF_1','-depsc')