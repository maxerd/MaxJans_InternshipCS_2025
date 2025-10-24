

figure(1001);clf
    bodemag(G_CL_dir);hold on;grid minor
    bodemag(G_CL_cp)
    bodemag(G_CL_2s)
    bodemag(G(1,1),'g--')
        xlim([1e-6 1e0])
        legend('CL, direct','CL, co-prime','CL, 2-stage','Model','location','south')

figure(1002);clf
    compare(sysID.data.CL.valid.direct,G_CL_dir,G_CL_cp,G_CL_2s,G(1,1));hold on;grid minor
        legend('Org. data','CL, direct','CL, co-prime','CL, 2-stage','Model','location','south')
