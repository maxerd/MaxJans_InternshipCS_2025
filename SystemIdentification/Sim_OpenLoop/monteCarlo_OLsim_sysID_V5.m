clear sys_frf_trd
clear sys_frf_lpm
clear sys_firstApprox
clear sys_initSys
clear sys_fixedOrder1
clear sys_fixedOrder5
clear sys_fixedOrder5
clear sys_freq

Nmc = 10;


compTrans = true;
newData   = true;
seperate_OLsim_sysID_V4

sys_frf_trd{1}     = sysID.nonPar.OL.frf.trd.filt;
sys_frf_lpm{1}     = sysID.nonPar.OL.frf.lpm.filt;
sys_frf_trdComp{1} = sysID.nonPar.OL.frfComp.trd.filt;
sys_frf_lpmComp{1} = sysID.nonPar.OL.frfComp.lpm.filt;
sys_firstApprox{1} = sysID.par.OL.time.firstApprox.filt;
sys_initSys{1}     = sysID.par.OL.time.initSys.prd.filt;
sys_fixedOrder1{1} = sysID.par.OL.time.fixedOrder1.prd.filt;
sys_fixedOrder5{1} = sysID.par.OL.time.fixedOrder5.prd.filt;
sys_straight{1}    = sysID.par.OL.time.straight.prd.filt;
sys_freq_lpm{1}    = sysID.par.OL.freq.lpm.prd.filt;
sys_freq_trd{1}    = sysID.par.OL.freq.trd.prd.filt;

for i=2:Nmc
    seperate_OLsim_sysID_V4

    sys_frf_trd{i}     = sysID.nonPar.OL.frf.trd.filt;
    sys_frf_lpm{i}     = sysID.nonPar.OL.frf.lpm.filt;
    sys_frf_trdComp{i} = sysID.nonPar.OL.frfComp.trd.filt;
    sys_frf_lpmComp{i} = sysID.nonPar.OL.frfComp.lpm.filt;
    sys_firstApprox{i} = sysID.par.OL.time.firstApprox.filt;
    sys_initSys{i}     = sysID.par.OL.time.initSys.prd.filt;
    sys_fixedOrder1{i} = sysID.par.OL.time.fixedOrder1.prd.filt;
    sys_fixedOrder5{i} = sysID.par.OL.time.fixedOrder5.prd.filt;
    sys_straight{i}    = sysID.par.OL.time.straight.prd.filt;
    sys_freq_lpm{i}    = sysID.par.OL.freq.lpm.prd.filt;
    sys_freq_trd{i}    = sysID.par.OL.freq.trd.prd.filt;
end

% save('_Data\monteCarlo_withTrans_8_Oct_25','sys_frf_trd','sys_frf_lpm','sys_firstApprox','sys_initSys','sys_fixedOrder1','sys_fixedOrder5','sys_straight','sys_freq_lpm','sys_freq_trd','sysID','G')
% save('_Data\monteCarlo_noTrans_8_Oct_25','sys_frf_trd','sys_frf_lpm','sys_firstApprox','sys_initSys','sys_fixedOrder1','sys_fixedOrder5','sys_straight','sys_freq_lpm','sys_freq_trd','sysID','G')

%%
figure(901);clf
for i=1:Nmc
    subplot(241)
        simpleBodemag(sys_firstApprox{i},'Hz',bodeRange);hold on
        simpleBodemag(G(1,1),'Hz','g--')
    subplot(242)
        simpleBodemag(sys_initSys{i},'Hz',bodeRange);hold on
        simpleBodemag(G(1,1),'Hz','g--')
    subplot(243)
        simpleBodemag(sys_fixedOrder1{i},'Hz',bodeRange);hold on
        simpleBodemag(G(1,1),'Hz','g--')
    subplot(244)
        simpleBodemag(sys_fixedOrder5{i},'Hz',bodeRange);hold on
        simpleBodemag(G(1,1),'Hz','g--')
    subplot(245)
        simpleBodemag(G(1,1)-sys_firstApprox{i},'Hz',bodeRange);hold on
    subplot(246)
        simpleBodemag(G(1,1)-sys_initSys{i},'Hz',bodeRange);hold on
    subplot(247)
        simpleBodemag(G(1,1)-sys_fixedOrder1{i},'Hz',bodeRange);hold on
    subplot(248)
        simpleBodemag(G(1,1)-sys_fixedOrder5{i},'Hz',bodeRange);hold on
end
subplot(221)
grid minor;xlim([bodeRange(1) bodeRange(end)])
subplot(222)
grid minor;xlim([bodeRange(1) bodeRange(end)])
subplot(223)
grid minor;xlim([bodeRange(1) bodeRange(end)])
subplot(224)
grid minor;xlim([bodeRange(1) bodeRange(end)])


figure(902);clf
for i=1:Nmc
    subplot(221)
        simpleBodemag(sys_frf_trdComp{i},'Hz');hold on
        simpleBodemag(G(1,1),'Hz','g--')
    subplot(222)
        simpleBodemag(sys_frf_lpmComp{i},'Hz');hold on
        simpleBodemag(G(1,1),'Hz','g--')
    subplot(223)
        simpleBodemag(G(1,1)-sys_frf_trdComp{i},'Hz');hold on
    subplot(224)
        simpleBodemag(G(1,1)-sys_frf_lpmComp{i},'Hz');hold on
end
subplot(211)
grid minor;xlim([bodeRange(1)*1000 bodeRange(end)])
subplot(212)
grid minor;xlim([bodeRange(1)*1000 bodeRange(end)])

figure(903);clf
for i=1:Nmc
    subplot(221)
        simpleBodemag(sys_freq_lpm{i},'Hz',bodeRange);hold on
        simpleBodemag(G(1,1),'Hz','g--',bodeRange)
    subplot(223)
        simpleBodemag(G(1,1)-sys_freq_lpm{i},'Hz',bodeRange);hold on
    subplot(222)
        simpleBodemag(sys_freq_trd{i},'Hz',bodeRange);hold on
        simpleBodemag(G(1,1),'Hz','g--',bodeRange)
    subplot(224)
        simpleBodemag(G(1,1)-sys_freq_trd{i},'Hz',bodeRange);hold on
end
subplot(211)
grid minor;xlim([bodeRange(1) bodeRange(end)])
subplot(212)
grid minor;xlim([bodeRange(1) bodeRange(end)])
























