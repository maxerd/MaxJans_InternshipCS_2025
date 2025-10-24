
figure(401);clf
sgtitle('Closed loop non-parametric identification of TMC setup')
    subplot(221)
        simpleBodemag(G_direct,'Hz');hold on;grid minor
        simpleBodemag(G_classical,'Hz')
        simpleBodemag(G_coprime_ry/G_coprime_ru,'Hz')
        simpleBodemag(G(1,1),'Hz','g--')
        xlim([1e-4 1e-1])
            title('Model estimations')
            legend('Direct method','Classical indirect method','Co-prime method','Original model','location','best')
    subplot(223)
        simpleBodephase(G_direct,'Hz','wrap');hold on;grid minor
        simpleBodephase(G_classical,'Hz','wrap')
        simpleBodephase(G_coprime_ry/G_coprime_ru,'Hz','wrap')
        simpleBodephase(G(1,1),'Hz','g--','wrap')
        xlim([1e-4 1e-1])
            title('Model estimations')
            legend('Direct method','Classical indirect method','Co-prime method','Original model','location','best')
    subplot(222)
        simpleBodemag(G_direct-G(1,1),'Hz');hold on;grid minor
        simpleBodemag(G_classical-G(1,1),'Hz')
        simpleBodemag(G_coprime_ry/G_coprime_ru-G(1,1),'Hz')
        xlim([1e-4 1e-1])
            title('Model estimation errors')
            legend('Direct method','Classical indirect method','Co-prime method','location','best')
    subplot(224)
        simpleBodephase(G_direct-G(1,1),'Hz','wrap');hold on;grid minor
        simpleBodephase(G_classical-G(1,1),'Hz','wrap')
        simpleBodephase(G_coprime_ry/G_coprime_ru-G(1,1),'Hz','wrap')
        xlim([1e-4 1e-1])
            title('Model estimation errors')
            legend('Direct method','Classical indirect method','Co-prime method','location','best')
