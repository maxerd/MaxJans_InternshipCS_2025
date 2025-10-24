

% Verification of G_{yr} using the residial test, only look at the cross correlation
figure(1613);
    resid(sysID.data.CL.valid.indirect1, est_yr)

% Verification of G_{ur} using the residial test, only look at the cross correlation
figure(1614);
    resid(sysID.data.CL.valid.indirect2, est_ur)

    
% Verify with previously obtainted model
figure(1615);
    bodemag(G_CL_cp,G(1,1));
    legend('Co-prime, closed loop','Non-param, open loop')
figure(1616);clf
    compare(sysID.data.CL.valid.direct,G_CL_cp)
figure(2616);clf
    resid(sysID.data.CL.valid.direct,G_CL_cp)






