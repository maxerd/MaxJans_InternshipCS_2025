
lumpedSystem_V3

%% Load in one of the implemented closed loop system measurements
meas = load('p__CB_none__ST_cl__SM_clw__SA_clw__DT_250903__MD_1h__WT_no__DS_1.mat');

%% Define the inputs
% tVec_C = 0:0.1:3600;
tVec_C = 0:0.1:(3600);

% Trajectory is the first 1/4 period of a cosine wave, followed by a
% constant. Giving a smooth ramp-like function
traj = 10.*[-0.5.*cos(2*pi/tVec_C(round(end/2))/2.*tVec_C(1:ceil(end/2)))+0.5 ones(1,10*tVec_C(ceil(end/2)))];

%% Simulate the system
y_CL = lsim(CL(:,1),traj,tVec_C);

%% Result visualization
startVal = 4030;
range = [startVal startVal+50];
% startVal = 0;
% range = [startVal startVal+3600];

figure(3001);clf
    subplot(211)
        plot(tVec_C,sqrt(y_CL(:,2)*190));hold on;grid minor
        % plot(meas.tVec+60,sqrt(meas.Watt*190),'m--')
        xlim(range)
            xlabel('Time [s]')
            ylabel('Heater voltage [V]')
            title('Voltage applied to heater')
            legend('Simulation','Control implementation','Location','best')
    subplot(212)
        plot(tVec_C,y_CL(:,2));hold on;grid minor
        % plot(meas.tVec+60,meas.Watt,'m--')
        % yline(0,'r--')
        xlim(range)
            xlabel('Time [s]')
            ylabel('Heater power [W]')
            title('Power send to heaters')
            legend('Simulated system','Implemented','Location','best')

shiftVal = 200;
errSig = traj(shiftVal:end)+22.1-meas.tempTM(1:(36000-shiftVal+2));
y_reverse = lsim(sysTMC([32 32+63 32+63+63],[3]),y_CL(:,2),tVec_C);
figure(3011);clf
    subplot(311)
        plot(tVec_C,traj+22.1,'r--');hold on;grid minor
        plot(tVec_C,y_CL(:,1)+22.1,'b')
        % plot(meas.tVec+shiftVal*meas.Ts,meas.tempTM,'g-.')
        xlim(range)
            xlabel('Time [s]')
            ylabel('Temperature [degC]')
            title('Temperature of Thermal Mass')
            legend('Trajectory','Simulated system','Implemented','Location','best')
    subplot(312)
        plot(tVec_C,traj-y_CL(:,1)');hold on;grid minor
        % plot(meas.tVec((shiftVal-1):36000)+shiftVal,traj(shiftVal:end)+22.1-meas.tempTM(1:(36000-shiftVal+2)),'m--')
        % plot(meas.tVec((shiftVal-1):36000)+shiftVal*meas.Ts,errSig,'m--')
        xlim(range)
            xlabel('Time [s]')
            ylabel('Tracking error [degC]')
            title('Tracking error of TM temperature')
            legend('Simulation error','Control implementation error','Location','best')
    subplot(313)
        plot(tVec_C,y_reverse+22.1);hold on;grid minor
        % plot(meas.tVec,meas.tempHT,'m--')
        xlim(range)
            xlabel('Time [s]')
            ylabel('Temperature [degC]')
            title('Temperature plots')
            legend('TM top simulation','TM bottom simulation','Heater simulation','Heater using control implementation','Location','best')


%% Error signal analysis
[PSD, CAS, freqVec] = fftCas_V2(errSig, 1/meas.Ts);

figure(3012);clf
    subplot(211)
        loglog(freqVec,PSD);grid minor
            xlabel('Frequency [Hz]')
            ylabel('Power [degC^2/Hz]')
            title('Power spectrum of tracking error')
    subplot(212)
        semilogx(freqVec,CAS);grid minor
            xlabel('Frequency [Hz]')
            ylabel('RMSe [degC]')
            title('CAS of tracking error')

%% Wattage signal analysis
[PSD, CAS, freqVec] = fftCas_V2(meas.Watt, 1/meas.Ts);
figure(3013);clf
    subplot(211)
        loglog(freqVec,PSD);grid minor
            xlabel('Frequency [Hz]')
            ylabel('Power [W^2/Hz]')
            title('Power spectrum of heater power')
    subplot(212)
        semilogx(freqVec,CAS);grid minor
            xlabel('Frequency [Hz]')
            ylabel('RMSe [W]')
            title('CAS of heater power')

%% Frequency domain plotting
opt = bodeoptions("cstprefs");
opt.phaseWrapping = 'off';
opt.FreqUnits = 'Hz';

SF = loopsens(d2c(contSys(:,1)),Ck.tot);
figure(3002);clf
    subplot(121)
        bode(SF.Si,SF.Ti,opt);grid minor
    subplot(122)
        bode(-SF.Li,opt);grid minor







