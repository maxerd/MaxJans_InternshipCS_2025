

% Verify with previously obtainted model
figure(617);clf;hold on;
    bode(G_CL_2s,G(1,1),opt)
    xlim([bodeRange(1) bodeRange(end)])
figure(1617);clf
    compare(sysID.data.CL.valid.direct,G_CL_2s)
figure(2617);clf
    resid(sysID.data.CL.valid.direct,G_CL_2s)