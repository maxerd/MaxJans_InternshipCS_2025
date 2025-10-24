

%% Unfiltered
ax7 = figure(baseFig_sim+7);clf
ax7.Name = 'comparison_bodes_withFRF';
set(gcf,'Position',[100 100 600 350])
    subplot(211);
        simpleBodemag(sysID.par.OL.freq.trd.prd.raw,'Hz',lw,bodeRange);hold on;grid minor
        simpleBodemag(sysID.par.OL.freq.lpm.prd.raw,'Hz',lw,bodeRange)
        simpleBodemag(sysID.nonPar.OL.frf.trd.raw,'Hz',1,'-.')
        simpleBodemag(sysID.nonPar.OL.frf.lpm.raw,'Hz',1,'-.')
        simpleBodemag(G(1,1),'Hz','g--',lw,bodeRange)
            title('Results from frequency domain identifications')
            % legend('nx=1','nx=2','nx=3','nx=4','nx=5','location','southwest')
            legend('Traditional','LPM','Traditional FRF','LPM FRF','Model','location','southwest')
    subplot(212);
        simpleBodemag(G(1,1)-sysID.par.OL.freq.trd.prd.raw,'Hz',lw,bodeRange);hold on;grid minor
        simpleBodemag(G(1,1)-sysID.par.OL.freq.lpm.prd.raw,'Hz',lw,bodeRange)
        simpleBodemag(G(1,1)-sysID.nonPar.OL.frf.trd.raw,'Hz',1,'-.')
        simpleBodemag(G(1,1)-sysID.nonPar.OL.frf.lpm.raw,'Hz',1,'-.')
            ylabel('Magnitude error [dB]')
            title('Magnitude error of frequency domain identifications')
            legend('Traditional','LPM','Traditional FRF','LPM FRF','location','southwest')
print('./_Figures/OL/parFreq/freq_comparison_bodes_withFRF','-depsc')

ax8 = figure(baseFig_sim+8);clf
ax8.Name = 'comparison_bodes_noFRF';
set(gcf,'Position',[100 100 600 350])
    subplot(211);
        simpleBodemag(sysID.par.OL.freq.trd.prd.raw,'Hz',lw,bodeRange);hold on;grid minor
        simpleBodemag(sysID.par.OL.freq.lpm.prd.raw,'Hz',lw,bodeRange)
        simpleBodemag(G(1,1),'Hz','g--',lw,bodeRange)
            title('Results from frequency domain identifications')
            % legend('nx=1','nx=2','nx=3','nx=4','nx=5','location','southwest')
            legend('Traditional','LPM','Model','location','southwest')
    subplot(212);
        simpleBodemag(G(1,1)-sysID.par.OL.freq.trd.prd.raw,'Hz',lw,bodeRange);hold on;grid minor
        simpleBodemag(G(1,1)-sysID.par.OL.freq.lpm.prd.raw,'Hz',lw,bodeRange)
            ylabel('Magnitude error [dB]')
            title('Magnitude error of frequency domain identifications')
            legend('Traditional','LPM','location','southwest')
print('./_Figures/OL/parFreq/freq_comparison_bodes_noFRF','-depsc')


%% Filtered
ax9 = figure(baseFig_sim+9);clf
ax9.Name = 'comparisonFilt_bodes_withFRF';
set(gcf,'Position',[100 100 600 350])
    subplot(211);
        simpleBodemag(sysID.par.OL.freq.trd.prd.raw,'Hz',lw,bodeRange);hold on;grid minor
        simpleBodemag(sysID.par.OL.freq.lpm.prd.raw,'Hz',lw,bodeRange)
        simpleBodemag(sysID.nonPar.OL.frf.trd.raw,'Hz',1,'-.')
        simpleBodemag(sysID.nonPar.OL.frf.lpm.raw,'Hz',1,'-.')
        simpleBodemag(G(1,1),'Hz','g--',lw,bodeRange)
            title('Results from frequency domain identifications')
            % legend('nx=1','nx=2','nx=3','nx=4','nx=5','location','southwest')
            legend('Traditional','LPM','Traditional FRF','LPM FRF','Model','location','southwest')
    subplot(212);
        simpleBodemag(G(1,1)-sysID.par.OL.freq.trd.prd.filt,'Hz',lw,bodeRange);hold on;grid minor
        simpleBodemag(G(1,1)-sysID.par.OL.freq.lpm.prd.filt,'Hz',lw,bodeRange)
        simpleBodemag(G(1,1)-sysID.nonPar.OL.frf.trd.filt,'Hz',1,'-.')
        simpleBodemag(G(1,1)-sysID.nonPar.OL.frf.lpm.filt,'Hz',1,'-.')
            ylabel('Magnitude error [dB]')
            title('Magnitude error of frequency domain identifications')
            legend('Traditional','LPM','Traditional FRF','LPM FRF','location','southwest')
print('./_Figures/OL/parFreq/freq_comparisonFilt_bodes_withFRF','-depsc')

ax10 = figure(baseFig_sim+10);clf
ax10.Name = 'comparisonFilt_bodes_noFRF';
set(gcf,'Position',[100 100 600 350])
    subplot(211);
        simpleBodemag(sysID.par.OL.freq.trd.prd.filt,'Hz',lw,bodeRange);hold on;grid minor
        simpleBodemag(sysID.par.OL.freq.lpm.prd.filt,'Hz',lw,bodeRange)
        simpleBodemag(G(1,1),'Hz','g--',lw,bodeRange)
            title('Results from frequency domain identifications')
            % legend('nx=1','nx=2','nx=3','nx=4','nx=5','location','southwest')
            legend('Traditional','LPM','Model','location','southwest')
    subplot(212);
        simpleBodemag(G(1,1)-sysID.par.OL.freq.trd.prd.filt,'Hz',lw,bodeRange);hold on;grid minor
        simpleBodemag(G(1,1)-sysID.par.OL.freq.lpm.prd.filt,'Hz',lw,bodeRange)
            ylabel('Magnitude error [dB]')
            title('Magnitude error of frequency domain identifications')
            legend('Traditional','LPM','location','southwest')
print('./_Figures/OL/parFreq/freq_comparisonFilt_bodes_noFRF','-depsc')

%% Fixed Orders unfiltered
ax11 = figure(baseFig_sim+11);clf
ax11.Name = 'freq_comparison_bodes_TRDfixedOrders';
% set(gcf,'Position',[100 100 600 350])
    subplot(211);
        simpleBodemag(sysID.par.OL.freq.trd1.prd.raw,'Hz',lw,bodeRange);hold on;grid minor
        simpleBodemag(sysID.par.OL.freq.trd2.prd.raw,'Hz',lw,bodeRange)
        simpleBodemag(sysID.par.OL.freq.trd3.prd.raw,'Hz',lw,bodeRange)
        simpleBodemag(sysID.par.OL.freq.trd4.prd.raw,'Hz',lw,bodeRange)
        simpleBodemag(sysID.par.OL.freq.trd5.prd.raw,'Hz',lw,bodeRange)
        simpleBodemag(sysID.nonPar.OL.frf.trd.raw,'Hz',1,'-.')
        simpleBodemag(G(1,1),'Hz','g--',lw,bodeRange)
            title('Results from different order frequency domain identifications, traditional')
            legend('nx=1','nx=2','nx=3','nx=4','nx=5','FRF','Model','location','southwest')
    subplot(212);
        simpleBodemag(G(1,1)-sysID.par.OL.freq.trd1.prd.raw,'Hz',lw,bodeRange);hold on;grid minor
        simpleBodemag(G(1,1)-sysID.par.OL.freq.trd2.prd.raw,'Hz',lw,bodeRange)
        simpleBodemag(G(1,1)-sysID.par.OL.freq.trd3.prd.raw,'Hz',lw,bodeRange)
        simpleBodemag(G(1,1)-sysID.par.OL.freq.trd4.prd.raw,'Hz',lw,bodeRange)
        simpleBodemag(G(1,1)-sysID.par.OL.freq.trd5.prd.raw,'Hz',lw,bodeRange)
        simpleBodemag(G(1,1)-sysID.nonPar.OL.frf.trd.raw,'Hz',1,'-.')
            ylabel('Magnitude error [dB]')
            title('Magnitude error of different order frequency domain identifications, traditional')
            legend('nx=1','nx=2','nx=3','nx=4','nx=5','FRF','location','southwest')
print('./_Figures/OL/parFreq/freq_comparison_bodes_TRDfixedOrders','-depsc')


ax12 = figure(baseFig_sim+12);clf
ax12.Name = 'freq_comparison_bodes_LPMfixedOrders';
% set(gcf,'Position',[100 100 600 350])
    subplot(211);
        simpleBodemag(sysID.par.OL.freq.lpm1.prd.raw,'Hz',lw,bodeRange);hold on;grid minor
        simpleBodemag(sysID.par.OL.freq.lpm2.prd.raw,'Hz',lw,bodeRange)
        simpleBodemag(sysID.par.OL.freq.lpm3.prd.raw,'Hz',lw,bodeRange)
        simpleBodemag(sysID.par.OL.freq.lpm4.prd.raw,'Hz',lw,bodeRange)
        simpleBodemag(sysID.par.OL.freq.lpm5.prd.raw,'Hz',lw,bodeRange)
        simpleBodemag(sysID.nonPar.OL.frf.lpm.raw,'Hz',1,'-.')
        simpleBodemag(G(1,1),'Hz','g--',lw,bodeRange)
            title('Results from different order frequency domain identifications, LPM')
            legend('nx=1','nx=2','nx=3','nx=4','nx=5','FRF','Model','location','southwest')
    subplot(212);
        simpleBodemag(G(1,1)-sysID.par.OL.freq.lpm1.prd.raw,'Hz',lw,bodeRange);hold on;grid minor
        simpleBodemag(G(1,1)-sysID.par.OL.freq.lpm2.prd.raw,'Hz',lw,bodeRange)
        simpleBodemag(G(1,1)-sysID.par.OL.freq.lpm3.prd.raw,'Hz',lw,bodeRange)
        simpleBodemag(G(1,1)-sysID.par.OL.freq.lpm4.prd.raw,'Hz',lw,bodeRange)
        simpleBodemag(G(1,1)-sysID.par.OL.freq.lpm5.prd.raw,'Hz',lw,bodeRange)
        simpleBodemag(G(1,1)-sysID.nonPar.OL.frf.lpm.raw,'Hz',1,'-.')
            ylabel('Magnitude error [dB]')
            title('Magnitude error of different order frequency domain identifications, LPM')
            legend('nx=1','nx=2','nx=3','nx=4','nx=5','location','southwest')
print('./_Figures/OL/parFreq/freq_comparison_bodes_LPMfixedOrders','-depsc')

%% Fixed Orders filtered
ax13 = figure(baseFig_sim+13);clf
ax13.Name = 'freq_comparison_bodes_TRDfixedOrders';
set(gcf,'Position',[100 100 600 350])
    subplot(211);
        simpleBodemag(sysID.par.OL.freq.trd1.prd.filt,'Hz',lw,bodeRange);hold on;grid minor
        simpleBodemag(sysID.par.OL.freq.trd2.prd.filt,'Hz',lw,bodeRange)
        simpleBodemag(sysID.par.OL.freq.trd3.prd.filt,'Hz',lw,bodeRange)
        simpleBodemag(sysID.par.OL.freq.trd4.prd.filt,'Hz',lw,bodeRange)
        simpleBodemag(sysID.par.OL.freq.trd5.prd.filt,'Hz',lw,bodeRange)
        simpleBodemag(sysID.nonPar.OL.frf.trd.filt,'Hz',1,'-.')
        simpleBodemag(G(1,1),'Hz','g--',lw,bodeRange)
            title('Results from different order frequency domain identifications, traditional')
            legend('nx=1','nx=2','nx=3','nx=4','nx=5','FRF','Model','location','southwest')
    subplot(212);
        simpleBodemag(G(1,1)-sysID.par.OL.freq.trd1.prd.filt,'Hz',lw,bodeRange);hold on;grid minor
        simpleBodemag(G(1,1)-sysID.par.OL.freq.trd2.prd.filt,'Hz',lw,bodeRange)
        simpleBodemag(G(1,1)-sysID.par.OL.freq.trd3.prd.filt,'Hz',lw,bodeRange)
        simpleBodemag(G(1,1)-sysID.par.OL.freq.trd4.prd.filt,'Hz',lw,bodeRange)
        simpleBodemag(G(1,1)-sysID.par.OL.freq.trd5.prd.filt,'Hz',lw,bodeRange)
        simpleBodemag(G(1,1)-sysID.nonPar.OL.frf.trd.filt,'Hz',1,'-.')
            ylabel('Magnitude error [dB]')
            title('Magnitude error of different order frequency domain identifications, traditional')
            legend('nx=1','nx=2','nx=3','nx=4','nx=5','FRF','location','southwest')
print('./_Figures/OL/parFreq/freq_comparisonFilt_bodes_TRDfixedOrders','-depsc')


ax14 = figure(baseFig_sim+14);clf
ax14.Name = 'freq_comparison_bodes_LPMfixedOrders';
set(gcf,'Position',[100 100 600 350])
    subplot(211);
        simpleBodemag(sysID.par.OL.freq.lpm1.prd.filt,'Hz',lw,bodeRange);hold on;grid minor
        simpleBodemag(sysID.par.OL.freq.lpm2.prd.filt,'Hz',lw,bodeRange)
        simpleBodemag(sysID.par.OL.freq.lpm3.prd.filt,'Hz',lw,bodeRange)
        simpleBodemag(sysID.par.OL.freq.lpm4.prd.filt,'Hz',lw,bodeRange)
        simpleBodemag(sysID.par.OL.freq.lpm5.prd.filt,'Hz',lw,bodeRange)
        simpleBodemag(sysID.nonPar.OL.frf.lpm.filt,'Hz',1,'-.')
        simpleBodemag(G(1,1),'Hz','g--',lw,bodeRange)
            title('Results from different order frequency domain identifications, LPM')
            legend('nx=1','nx=2','nx=3','nx=4','nx=5','FRF','Model','location','southwest')
    subplot(212);
        simpleBodemag(G(1,1)-sysID.par.OL.freq.lpm1.prd.filt,'Hz',lw,bodeRange);hold on;grid minor
        simpleBodemag(G(1,1)-sysID.par.OL.freq.lpm2.prd.filt,'Hz',lw,bodeRange)
        simpleBodemag(G(1,1)-sysID.par.OL.freq.lpm3.prd.filt,'Hz',lw,bodeRange)
        simpleBodemag(G(1,1)-sysID.par.OL.freq.lpm4.prd.filt,'Hz',lw,bodeRange)
        simpleBodemag(G(1,1)-sysID.par.OL.freq.lpm5.prd.filt,'Hz',lw,bodeRange)
        simpleBodemag(G(1,1)-sysID.nonPar.OL.frf.lpm.filt,'Hz',1,'-.')
            ylabel('Magnitude error [dB]')
            title('Magnitude error of different order frequency domain identifications, LPM')
            legend('nx=1','nx=2','nx=3','nx=4','nx=5','location','southwest')
print('./_Figures/OL/parFreq/freq_comparisonFilt_bodes_LPMfixedOrders','-depsc')

%%

preds_freqApproaches = compare(sysID.OLdataFilt.valid.id,...
                            sysID.par.OL.freq.trd.prd.filt,...
                            sysID.par.OL.freq.lpm.prd.filt);
sims_freqApproaches  = compare(sysID.OLdataFilt.valid.id,...
                           sysID.par.OL.freq.trd.prd.filt,...
                           sysID.par.OL.freq.lpm.prd.filt);

preds_freqApproaches_err{1} = sysID.OLdataFilt.valid.y-squeeze(preds_freqApproaches{1}.y);
preds_freqApproaches_err{2} = sysID.OLdataFilt.valid.y-squeeze(preds_freqApproaches{2}.y);

sims_freqApproaches_err{1} = sysID.OLdataFilt.valid.y-squeeze(sims_freqApproaches{1}.y);
sims_freqApproaches_err{2} = sysID.OLdataFilt.valid.y-squeeze(sims_freqApproaches{2}.y);

disp(' ')
disp('RMS prediction errors of time based identification approaches:')
disp(['LPM based, RMSE: ',  num2str(rms(preds_freqApproaches_err{1}),3),'degC'])
disp(['Traditional, RMSE: ',num2str(rms(preds_freqApproaches_err{2}),3),'degC'])
disp(' ')
disp('RMS simulation errors of time based identification approaches:')
disp(['LPM based, RMSE: ',  num2str(rms(sims_freqApproaches_err{1}),3),'degC'])
disp(['Traditional, RMSE: ',num2str(rms(sims_freqApproaches_err{2}),3),'degC'])


% ax15 = figure(baseFig_sim+15);clf
% ax15.Name = 'freq_comparison_compare_LPMfixedOrders';
% set(gcf,'Position',[100 100 600 350])
%     plot(sysID.OLdata.valid.tVec,sysID.OLdataFilt.valid.y);hold on
%     plot(sysID.OLdata.valid.tVec,preds_freqApproaches{1}.y);hold on
%     plot(sysID.OLdata.valid.tVec,preds_freqApproaches{2}.y)
















