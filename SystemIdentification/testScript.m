
trainData = sysID.OLdata.train;
timeData = sysID.OLdata.train.tVec;
% trainData = sysID.CLdata.train.direct;
% timeData = sysID.CLdata.train.tVec;


tic
simOut = lsim(sysID.par.OL.time.firstApprox.filt,trainData.u,timeData);
% simOut = lsim(sysID.par.OL.time.initSys.sim.filt,trainData.u,timeData);
% diff = iddata(trainData.y-simOut,trainData.u);

detrendOut = fit(timeData',trainData.y,'exp2');
func = @(x) detrendOut.a*exp(detrendOut.b*x)+detrendOut.c*exp(detrendOut.d*x);
newy = func(timeData)';
diff = iddata(trainData.y-newy+mean(trainData.y(end-1000:end)),trainData.u,Ts);


[testFRF, ~] = makeOpenLoopFRF_sysIdent(trainData.y-newy+mean(trainData.y(end-1000:end)), trainData.u, fs);


figure(1);clf;plot(diff)

newSys = ssest(diff);%+sysID.par.OL.time.firstApprox.filt;
t1 = toc;
tic
ssest(sysID.OLdata.train.id,init_sys);
t2 = toc;

% trainData = sysID.CLdata.train.direct;
% timeData = sysID.CLdata.train.tVec;

% simOut = lsim(newSys,trainData.u,timeData);
% diff = iddata(trainData.y-simOut,trainData.u);
% 
% newSys = ssest(diff,1)+newSys;

%%
figure(2);clf
    subplot(211)
        bodemag(newSys);hold on
        % bodemag(sysID.par.OL.time.firstApprox.raw)
        bodemag(testFRF)
        bodemag(sysID.nonPar.OL.frf.trd.raw)
        bodemag(sysID.nonPar.OL.frf.lpm.raw)
        bodemag(G(1,1),bodeRange,opt)
    subplot(212)
        bodemag(G(1,1)-newSys);hold on
        bodemag(G(1,1)-sysID.par.OL.time.firstApprox.raw)

%%
figure(3);clf
    % compare(sysID.OLdata.train.id,sysID.par.OL.time.firstApprox.filt+newSys,sysID.par.OL.time.firstApprox.filt,G(1,1));hold on
    % compare(iddata(sysID.CLdataFilt.valid.y,sysID.CLdata.valid.u,10),newSys,sysID.par.OL.time.firstApprox.filt,sysID.par.OL.time.fixedOrder5.prd.filt);hold on
    compare(sysID.OLdataFilt.valid.id,...
            newSys,...
            sysID.par.OL.time.firstApprox.raw,...
            sysID.par.OL.time.fixedOrder1.prd.raw,...
            sysID.par.OL.time.fixedOrder5.prd.raw,...
            sysID.par.OL.time.initSys.prd.raw,...
            sysID.par.OL.time.straight.prd.raw);hold on

%%
figure(4);clf
    % compare(sysID.OLdata.train.id,sysID.par.OL.time.firstApprox.filt+newSys,sysID.par.OL.time.firstApprox.filt,G(1,1));hold on
    resid(iddata(sysID.CLdata.full.y,sysID.CLdata.full.u,10),newSys,sysID.par.OL.time.firstApprox.filt,sysID.par.OL.time.fixedOrder5.prd.filt);hold on
    % resid(sysID.OLdata.valid.id,...
    %         newSys,...
    %         sysID.par.OL.time.firstApprox.raw,...
    %         sysID.par.OL.time.initSys.prd.raw);hold on



