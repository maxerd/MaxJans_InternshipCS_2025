
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
    % print('./_Figures/OL/parTime/firstOrderApprox_1','-depsc')
