
[sysID.OLdata, sysID.OLdataFilt] = makeOLdata(G,tMeas,dist,Ts);

trainData = sysID.OLdata.train;
timeData = sysID.OLdata.train.tVec;

DCgain = mean(trainData.y(end-1000:end));

%% Using first order approx
firstApproxSim      = lsim(sysID.par.OL.time.firstApprox.filt,trainData.u,timeData);
firstApproxSim_diff = iddata(trainData.y-firstApproxSim+DCgain,trainData.u,Ts);
figure(1);plot(firstApproxSim_diff)
figure(11);clf;hold on
    plot(timeData,trainData.y)
    plot(timeData,firstApproxSim)

[testFRF_first, ~] = makeOpenLoopFRF_sysIdent(firstApproxSim_diff.y+DCgain, trainData.u, fs);

firstSys = ssest(firstApproxSim_diff,1);

figure(9);bode(firstSys,testFRF_first,G(1,1),'g--',opt);grid minor

%% Using exponential approx

[sysID.OLdata, sysID.OLdataFilt] = makeOLdata(G,tMeas,dist,Ts);

trainData = sysID.OLdata.train;
timeData = sysID.OLdata.train.tVec;

DCgain = mean(trainData.y(end-1000:end));

detrendOut = fit(timeData',trainData.y,'exp2');
func = @(x) detrendOut.a*exp(detrendOut.b*x)+detrendOut.c*exp(detrendOut.d*x);
newy = func(timeData)';
expApproxSim_diff = iddata(trainData.y-newy+mean(trainData.y(end-1000:end)),trainData.u,Ts);
figure(2);plot(expApproxSim_diff)

[testFRF_exp, ~] = makeOpenLoopFRF_sysIdent(expApproxSim_diff.y+DCgain, trainData.u, fs);

expSys = ssest(expApproxSim_diff);

figure(10);
    subplot(211)
        bodemag(expSys,...
                sysID.par.OL.time.firstApprox.raw,...
                testFRF_exp,...
                sysID.par.OL.time.fixedOrder5.prd.raw,...
                G(1,1),'g--',opt);grid minor
            legend('Par, exp removed','First order approx','FRF','Model')
    subplot(212)
        bodemag(G(1,1)-expSys,...
                G(1,1)-sysID.par.OL.time.firstApprox.raw,...
                G(1,1)-sysID.par.OL.time.fixedOrder5.prd.raw,...
                G(1,1)-testFRF_exp,opt);grid minor
            legend('Par, exp removed','First order approx','FRF')

figure(110);
    compare(sysID.OLdata.valid.id,expSys,sysID.par.OL.time.firstApprox.raw,sysID.par.OL.time.initSys.prd.raw);grid minor



