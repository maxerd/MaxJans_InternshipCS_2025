
%% Calculate and plot the first order approximation
steady_val = mean(sysID.OLdata.full.y(end-100:end));
tau = find(sysID.OLdata.full.y>0.632.*steady_val(end),1,'first');

figure(baseFig_sim+1);clf
set(gcf,'position',[100 100 800 350])
    subplot(221);hold on
        plot(sysID.OLdata.full.tVec,sysID.OLdata.full.y);grid minor
        xline(tau*Ts,'g--',lineWidth=lw)
        yline(steady_val,'r--',lineWidth=lw)
        yline(sysID.OLdata.full.y(tau),'m--',lineWidth=lw)
        plot(tau*Ts,sysID.OLdata.full.y(tau),'ro',lineWidth=2)
            xlabel('Time [s]')
            ylabel('Temperature [degC]')
            title('Step response of thermal system')
            legend('Temperature',['Time constant, $\tau$: ',num2str(tau*Ts),'s'],['Steady state value, $A$: ',num2str(steady_val),'degC'],['Time constant cut: ',num2str(sysID.OLdata.full.y(tau)),'degC'],'Location','southeast')
    subplot(223);hold on
        plot(sysID.OLdata.full.tVec,sysID.OLdata.full.u);grid minor
        yline(mean(sysID.OLdata.full.u),'r--',LineWidth=lw)
        xlim([sysID.OLdata.full.tVec(1) sysID.OLdata.full.tVec(end)])
            xlabel('Time [s]')
            ylabel('Input power [W]')
            title('Input power used for step response identification')
            legend('Input signal',['Mean wattage: ',num2str(mean(sysID.OLdata.full.u)),'W'])
    subplot(222)
        simpleBodemag(sysID.par.OL.time.firstApprox.raw,'Hz',lw,bodeRange);hold on;grid minor
        simpleBodemag(G(1,1),'Hz','k--',lw,bodeRange)
        xline(1/tau/Ts,'g--',lineWidth=lw)
        yline(mag2db(steady_val/mean(sysID.OLdata.full.u)),'r--',lineWidth=lw)
        xlim([bodeRange(1) bodeRange(end)])
            title('Bode magnitude identified from step response')
            legend('Identified system','Model',['Cutt-off frequency, $1/\tau$: ',num2str(1/tau/Ts*10000,3),'$\cdot 10^{-4}$ Hz'],['DC-gain, $\frac{A}{W_m}$: ',num2str(mag2db(evalfr(sysID.par.OL.time.firstApprox.raw,0)),3),'dB'],'Location','southwest')
    subplot(224)
        simpleBodephase(sysID.par.OL.time.firstApprox.raw,'Hz',lw,bodeRange,'wrap');hold on;grid minor
        simpleBodephase(G(1,1),'Hz','k--',lw,bodeRange,'wrap')
        xline(1/tau/Ts,'g--',lineWidth=lw)
        xlim([bodeRange(1) bodeRange(end)])
print('./_Figures/OL/parTime/firstOrderApprox','-depsc')