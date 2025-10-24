
figure(611);clf;hold on;
sgtitle('Bode plots from direct methods');
    subplot(221)
        simpleBodemag(G_CL_dir,'Hz',lw);hold on;grid minor
        simpleBodemag(G(1,1),'Hz',lw);
            xlim([bodeRange(1) bodeRange(end)])
            xline(fs/2,LineWidth=1)
            title('Plant model');
            legend('Closed loop, direct method','Model','Location','best');
    subplot(223)
        simpleBodephase(G_CL_dir,'Hz',lw,'wrap');hold on;grid minor
        simpleBodephase(G(1,1),'Hz',lw,'wrap');
            xlim([bodeRange(1) bodeRange(end)])
            xline(fs/2,LineWidth=1)
            title('Plant model');
            legend('Closed loop, direct method','Model','Location','best');
    subplot(222)
        simpleBodemag(H_CL_dir,'Hz',lw);hold on;grid minor
        % simpleBodemag(H,'Hz',lw)
            xline(fs/2,LineWidth=1)
            title('Plant model');
            legend('Closed loop, direct method','Location','best');
    subplot(224)
        simpleBodephase(H_CL_dir,'Hz',lw,'wrap');hold on;grid minor
            xline(fs/2,LineWidth=1)
            title('Plant model');
            legend('Closed loop, direct method','Location','best');

figure(1611);clf
    compare(sysID.data.CL.valid.direct,G_CL_dir)
figure(1612);clf
    resid(sysID.data.CL.valid.direct,G_CL_dir)