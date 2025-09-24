clc

addpath(genpath("..\Consts\"))
addpath(genpath(".\Funcs\"))

%% User definitions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Options for the user are shown in this section %%%%%%%%%


% Make and save figures?
    % -1 -> All figures, 
    %  0 -> No figures, 
    %  1 -> Static lumps visualization, 
    %  2 -> Bode plots, 
    %  3 -> Simulations, 
    %  4 -> Controller simulations, 
    %  5 -> Dynamic lumps visualization
    % makeFigs = [3]; 
    % makeFigs = [3 5]; 
    % makeFigs = [5]; 
    makeFigs = [0]; 
    saveFigs = false;

% Define what the system output will be
    % -1-> All states, 
    %  0-> Only top middle lump, 
    %  1-> 0+Heater temp
    outputOpt = -1; 

% Define what the system input will be
    % -1 -> All heaters, 
    %  0 -> None of the heaters, 
    %  1 -> Only left most heater (HT1), 
    %  2 -> Only second to most left heater (HT2), 
    %  3 -> Only middle heater (HT3), 
    %  4 -> Only second to most right heater (HT4), 
    %  5 -> Only right most heater (HT5)
    inputOpt = 3; 
    inputAmp = 40; % Input amplitude (Watt)

% Make and simulate controller?
    makeController = false;
    bwController   = 0.0025;
    saveController = false;

% Define the (constant) water flow
    % waterVelocity = 0; %[ml/s]
    waterVelocity = ureal('waterVel',0.5,'PlusMinus',[-0.5 0.5]); %[ml/s]

% Make LFR repesentation
    makeLFR = true; %[true/false]

%%%%%%%% Normally, no changes past this point are needed %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Define constants
% Load in some thermal and mechanical constants
c = constsGen('thermal','mech');

% The TM is divided into 9 lumps in the x direction
%                        7 lumps in the z direction
%                        3 lumps in the y direction
% This gives a total of 63 lumps per layer (layer = per y position)
Nx = 9;
Nz = 7;
Ny = 3;
Nlumps = Nx*Nz;

% A lot of data is aquired on the lumps, and manually put into a file for
% matlab acces, make sure the version corresponds to the current script version
manualLumpData_water_V3

%% Make graph L, A and D matrices
% Using graph theory (ModRed_6a slides, Model Reduction course TU/e), take
% the lumps and connect them using the conductances as defined in 'manualLumpData_water_V3'
makeGraphMatrices_water

%% Make state space E matrix
% Make the E matrix (Determine the thermal masses)
makeEmatrix_water

%% Make state space C matrix
% Make the C matrix (Determine the system outputs)
switch outputOpt
    case -1
        C = eye(Ny*Nlumps+29);
    case 0
        C = zeros(1,Ny*Nlumps+29);C(1,32) = 1;
    case 1
        C = zeros(3,Ny*Nlumps+29);C(1,32) = 1;C(2,32+63) = 1;C(3,end) = 1;
end

try
    C;
catch
    error('Invalid outputOpt')
end

%% System definition and basic analysis
% clc

% Make a descriptor state space model
    % This can be done using a dss system (descriptor state space), however,
    % this can give problems if the thermal masses are incorrectly defined.
    % The system can also be written as a 'standard' state space system,
    % using the inverse of E

    % sysTMC = dss(-L,B_heaters,C,0,E);
    sysTMC = ss(inv(E)*-L,inv(E)*B_heaters,C,0);
    % sysTMC_w = ss(inv(E)*-L,inv(E)*B_heaters,C,0);

% Some basic analysis of the system (bode plots, obsv, ctrb)
    % This script can be expanded in the future
    systemAnalysis

%% System similation initialization
% Define the desired time vector for the simulation
    % tVec = 0:10:(4*3600);
    tVec = 0:10:(8*3600);
    % inputAmp = 40.*(sin(0.001*tVec)+1);

% Different input signal options
    switch inputOpt
        case -1
            u = chooseInput(tVec,'all',inputAmp);
        case 0
            u = chooseInput(tVec,'none',0);
        case 1
            u = chooseInput(tVec,'HT1',inputAmp);
        case 2
            u = chooseInput(tVec,'HT2',inputAmp);
        case 3
            u = chooseInput(tVec,'HT3',inputAmp)';
        case 4
            u = chooseInput(tVec,'HT4',inputAmp);
        case 5
            u = chooseInput(tVec,'HT5',inputAmp);

        case 6
            u = chooseInput(tVec,'HT1',inputAmp,'HT3',inputAmp,'HT5',inputAmp);
    end
    if size(u,1)==5
        u = u';
    end


% Different frequently used initial condition options
    x0 = 21.7.*ones(Ny*Nlumps+29,1);
    % x0 = 20.6.*ones(Ny*Nlumps+1,1);
    % x0 = 20.*ones(Ny*Nlumps+1,1);
    % x0 = 94.5.*ones(Ny*Nlumps+29,1);

% Different frequently used ambient temperature options
    Tamb = 21.7.*ones(length(u),1);
    % Tamb = 20.6.*ones(length(u),1);
    % Tamb = 20.*ones(length(u),1);

% Different frequently used initial water temperature options
    % Twat = 21.7.*ones(length(u),1);
    % Twat = 20.6.*ones(length(u),1);
    Twat = 20.*ones(length(u),1);

%% Linear system simulation
% Using the previously defined variables

clear tempSimData xSimData

[tempSimData,~,xSimData]  = lsim(sysTMC,[u Tamb Twat],tVec,x0);

%% Simulation visualization

if or(makeFigs==-1,nonzeros(makeFigs==3))
    
% plot of the heater lumps
    figure(3001);clf
        sgtitle('Temperature of heater lumps')
    for n=0:8
        subplot(331+n)
            plot(tVec,tempSimData(:,[(127):(127+6)]+n*7),LineWidth=1);grid minor
                xlabel('Time [s]')
                ylabel('Temperature [degC]')
                title(['Heater lump ',num2str(127+n*7),' to ',num2str(127+6+n*7)])
    end

% Plot of TM temperatures 
    figure(3002);clf;hold on;grid minor
        plot(tVec,tempSimData(:,[32 32+63 32+63+63]),LineWidth=1)
            legend('TM top simulation','TM bottom simulation','Heater simulation','Location','best')
            xlabel('Time [s]')
            ylabel('Temperature [degC]')
            title('Temperature of selected Thermal mass lumps')
   
% Plot of TM temperatures including water
    figure(3003);clf;hold on;grid minor
        plot(tVec,tempSimData(:,[32 32+63 32+63+63 63+63+63+1 63+63+63+29]),LineWidth=1)
        legend('TM top simulation','TM bottom simulation','Heater simulation','Ingoing water simulation','Outgoing water simulation','Location','best')
            xlabel('Time [s]')
            ylabel('Temperature [degC]')
            title('Temperature of Thermal mass')

% Plot of all water lumps
    figure(3004);clf;hold on;grid minor
        plot(tVec,tempSimData(:,63+63+63+[1:29]),LineWidth=1)
            xlabel('Time [s]')
            ylabel('Temperature [degC]')
            title('Temperature of water lumps')

% Plot of all TM top lumps
    figure(3005);clf;hold on;grid minor
        plot(tVec,tempSimData(:,[1:126]),LineWidth=1)
            xlabel('Time [s]')
            ylabel('Temperature [degC]')
            title('Temperature of Thermal mass lumps')

% Calculate the length of the water path
    waterLengths = zeros(1,15);
    waterLengths(1)     = LTM6.Lx/2+5;
    waterLengths(2)     = LTM13.Lx;
    waterLengths(3)     = LTM20.Lx;
    waterLengths(4:6)   = LTM27.Lx;
    waterLengths(7:8)   = (LTM48.Lx/2+5)+(LTM48.Ly/2+5);
    waterLengths(9:11)  = LTM27.Lx;
    waterLengths(12:13) = (LTM48.Lx/2+5)+(LTM48.Ly/2+5);
    waterLengths(14:15) = LTM27.Lx;
    waterLengths = [waterLengths fliplr(waterLengths(1:end-1))]*triu(ones(29));

% Water lump temperature
    figure(3006);clf;hold on;grid minor
        plot(waterLengths/(waterVelocity*1e9/(pi*5^2)),tempSimData(end,63+63+63+[1:29]),'x-');
            xlabel('Time in TM [sec]');
            ylabel('Temperature [degC]');
            title(['Water temperature over time, flow: ',num2str(waterVelocity*1e6),'ml/s'])

% Additional water lump temperature plots, with respect to lump and distance
    % figure(3007);clf;hold on;grid minor
    %     plot(tempSimData(end,63+63+63+[1:29]),'x-');
    %         xlabel('Water lump [-]');
    %         ylabel('Temperature [degC]');
    %         title('Water temperature per lump')
    % figure(3008);clf;hold on;grid minor
    %     plot(waterLengths,tempSimData(end,63+63+63+[1:29]),'x-');
    %         xlabel('Distance travelled in TM [mm]');
    %         ylabel('Temperature [degC]');
    %         title('Water temperature over travelled distance')


end

%% Sanity checks
% A few lines that can be used to manually see if the model makes sense
% Normally running this script is unnessecary

% sanityChecks

%% Make LFR repesentation
inputNames = {'HT1'};
for i=2:(length(sysTMC.inputName)-2)
    inputNames = [inputNames, (['HT' num2str(i)])];
end
inputNames = [inputNames, 'tempAmb', 'tempWater_in'];
sysTMC.inputName = inputNames;

outputNames = {'Tmp1'};
for i=2:(length(sysTMC.outputName))
    outputNames = [outputNames, (['Tmp' num2str(i)])];
end
sysTMC.outputName = outputNames;

if makeLFR
    G = sysTMC;

    refSum_P = sumblk('tempErr = tempRef - Tmp32');

    P = connect(G,refSum_P,{'tempRef','tempAmb','tempWater_in',inputNames{1:5}},{outputNames{1:end},inputNames{1:5},'tempErr'});
    Pnu = 5;
    Pny = 1;

    % An Hinf controller or something similar can be defined here using the
    % LFR representation and the Pnu (amount of inputs) and Pny (amount of outputs)
    if strcmp(class(waterVelocity),'umat')
        [Pdel,Delta,Blkstruct] = lftdata(P);
        Pdel.outputName{1} = 'z_waterVel_delta';
        Pdel.inputName{1}  = 'w_waterVel_delta';

        Dnu = 1;
        Dny = 1;
    end
else % Controller generation only works for non-uncertain system

%% Generate simple PID controller for SISO system (P_heater --> Top middle lump)

% Select the system that needs to be controlled
% contSys = newSys;
contSys = c2d(sysTMC(32,[3 6]),0.001);

% Define the system I/O
contSys.InputName =  {'W','Tamb'};
contSys.OutputName = 'T';

if makeController
% Generate a PID controller for a user selected bandwidth
    Ck = genPID(d2c(contSys(1,1)),bwController);
    Ck.tot.InputName =  'Terr';
    Ck.tot.OutputName = 'W';
    
    refSum = sumblk('Terr = Tref - T');
    
% Make the closed loop connection
    CL = connect(d2c(contSys),Ck.tot,refSum,{'Tref','Tamb'},{'T','W'});
    
% Save the controller
    if saveController
        save(['../Controller/PID_',datestr(now, 'yymmdd_HHMMSS'),'.mat'],'Ck','contSys','CL')
    end

%% Simulate the closed loop system (including the generated controller)

    if or(makeFigs==-1,nonzeros(makeFigs==4))
        
    % Define the inputs
        tVec_C = 0:0.1:3600;
    
    % Trajectory is the first 1/4 period of a cosine wave, followed by a
    % constant. Giving a smooth ramp-like function
        traj = 10.*[-0.5.*cos(2*pi/tVec_C(round(end/2))/2.*tVec_C(1:ceil(end/2)))+0.5 ones(1,10*tVec_C(ceil(end/2)))];
    
    % Simulate the system
        y_CL = lsim(CL(:,1),traj,tVec_C);

    % Save the controller
        if saveController
            save(['../Controller/PID_',datestr(now, 'yymmdd_HHMMSS'),'.mat'],'Ck','contSys','CL','traj','tVec_C')
        end
        
        figure(4001);clf
            subplot(211)
                plot(tVec_C,sqrt(y_CL(:,2)*190));hold on;grid minor
                    xlabel('Time [s]')
                    ylabel('Heater voltage [V]')
                    title('Voltage applied to heater')
                    legend('Simulation','Location','best')
            subplot(212)
                plot(tVec_C,y_CL(:,2));hold on;grid minor
                yline(0,'r--')
                    xlabel('Time [s]')
                    ylabel('Heater power [W]')
                    title('Power send to heaters')
                    legend('Simulated system','Minimum power','Location','best')
        
        y_reverse = lsim(sysTMC([32 32+63 32+63+63],[3]),y_CL(:,2),tVec_C);
        figure(4002);clf
            subplot(311)
                plot(tVec_C,traj+22.1,'r--');hold on;grid minor
                plot(tVec_C,y_CL(:,1)+22.1,'b')
                    xlabel('Time [s]')
                    ylabel('Temperature [degC]')
                    title('Temperature of Thermal Mass')
                    legend('Trajectory','Simulated system','Location','best')
            subplot(312)
                plot(tVec_C,traj-y_CL(:,1)');hold on;grid minor
                    xlabel('Time [s]')
                    ylabel('Tracking error [degC]')
                    title('Tracking error of TM temperature')
                    legend('Simulation error','Location','best')
            subplot(313)
                plot(tVec_C,y_reverse+22.1);hold on;grid minor
                    xlabel('Time [s]')
                    ylabel('Temperature [degC]')
                    title('Temperature plots')
                    legend('TM top simulation','TM bottom simulation','Heater simulation','Location','best')
    
    end

end

end
%% Dyamic lump visualization

if or(makeFigs==-1,nonzeros(makeFigs==5))
    if ~nonzeros(makeFigs==3)
        warning('No new simulation made, using old data')
    end
    try
        TM_Visual(tempSimData,tVec)
    catch
        error('Problem in dynamic lump visualization')
    end
end




