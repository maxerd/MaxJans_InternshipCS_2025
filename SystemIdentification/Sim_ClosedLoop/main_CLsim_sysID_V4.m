clear all
% close all
clc

%% Pre-requisites
opt = bodeoptions('cstprefs');
opt.PhaseWrapping = 'on';
opt.FreqUnits = 'Hz';

lw = 1.5; % Define the wanted linewidth for the plots
bodeRange = logspace(-6,-1,100);

list_factory = fieldnames(get(groot,'factory'));
index_interpreter = find(contains(list_factory,'Interpreter'));
for i = 1:length(index_interpreter)
    default_name = strrep(list_factory{index_interpreter(i)},'factory','default');
    set(groot, default_name,'latex');
end

%% Main variables

% Figure number of of which all others are based
baseFig_sim = 8000;

% Save the identification?
saveID = true; % [true/false]

% Total simulation time
tMeas = 20; % [h]
% Sampling time of simulation
Ts = 10; % [h]
fs = 1/Ts; % [h]

%% Make the system used for data generation
disp('Constructing system model')
main_lumpedSystem_water_V3

G = sysTMC(32,[3 6]);

%% Define some controller variables
controllerBW = 0.0001; % [Hz] Controller bandwidth

%% Define the identification signal
% maxAmplitude = 24; % [W] Maximum exitation signal
maxAmplitude = 120; % [W] Maximum exitation signal
positiveOnly = 1;  % [-] Define whether the exitation signal can be negative

disp('Generating multisine for identification')
dist     = genMultisine(fs, tMeas*3600*fs, 1, maxAmplitude,     positiveOnly)+25;

%% Run simulation and format data for later use
% [sysID.data, sysID.dataFilt]= makeCLdata(G,tMeas,controllerBW,dist,Ts);
% [~, sysID.data] = makeCLdata(G,tMeas,controllerBW,dist,Ts);
[~, sysID.data,C] = makeCLdata(G,tMeas,controllerBW,dist,Ts);

%% Do a non-parametric identification
[G_direct,~]  = makeClosedLoopFRF(sysID.data.CL.train.e,sysID.data.CL.train.r,C,fs);
[G_classical,~] = makeOpenLoopFRF_sysIdent(sysID.data.CL.train.y,sysID.data.CL.train.u,fs);
[G_coprime_ru,~] = makeOpenLoopFRF_sysIdent(sysID.data.CL.train.u,sysID.data.CL.train.r,fs);
[G_coprime_ry,~] = makeOpenLoopFRF_sysIdent(sysID.data.CL.train.y,sysID.data.CL.train.r,fs);

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
    subplot(222)
        simpleBodemag(G_direct-G(1,1),'Hz');hold on;grid minor
        simpleBodemag(G_classical-G(1,1),'Hz')
        simpleBodemag(G_coprime_ry/G_coprime_ru-G(1,1),'Hz')
        xlim([1e-4 1e-1])
            title('Model estimation errors')
            legend('Direct method','Classical indirect method','Co-prime method','location','best')

%% Do a parametric identification, direct method
sysID.data.CL.train.direct = iddata(sysID.data.CL.train.y,sysID.data.CL.train.u,Ts);
sysID.data.CL.valid.direct = iddata(sysID.data.CL.valid.y,sysID.data.CL.valid.u,Ts);
sysID.data.CL.trans.direct = iddata(sysID.data.CL.trans.y,sysID.data.CL.trans.u,Ts);

% structure
nb = 2;
nc = 1;
nd = 1;
nf = 3;
nk = 0;
est = bj(sysID.data.CL.train.direct, [nb nc nd nf nk]);

% Make the transfer function with the identified parameters
G_CL_dir = tf(est.B,est.F,est.Ts,'variable','q^-1');
H_CL_dir = tf(est.C,est.D,est.Ts,'variable','q^-1');

figure(611);clf;hold on;
sgtitle('Bode plots from direct methods');
    subplot(221)
        simpleBodemag(G_CL_dir,'Hz',lw);hold on;grid minor
        simpleBodemag(G(1,1),'Hz',lw);
            xlim([bodeRange(1) bodeRange(end)])
            xline(fs/2,LineWidth=1)
            title('Plant model');
            legend('Closed loop, direct method','Model','Location','best');
    subplot(223)
        simpleBodephase(G_CL_dir,'Hz',lw,'wrap');hold on;grid minor
        simpleBodephase(G(1,1),'Hz',lw,'wrap');
            xlim([bodeRange(1) bodeRange(end)])
            xline(fs/2,LineWidth=1)
            title('Plant model');
            legend('Closed loop, direct method','Model','Location','best');
    subplot(222)
        simpleBodemag(H_CL_dir,'Hz',lw);hold on;grid minor
        simpleBodemag(H,'Hz',lw)
            xline(fs/2,LineWidth=1)
            title('Plant model');
            legend('Closed loop, direct method','Location','best');
    subplot(224)
        simpleBodephase(H_CL_dir,'Hz',lw,'wrap');hold on;grid minor
            xline(fs/2,LineWidth=1)
            title('Plant model');
            legend('Closed loop, direct method','Location','best');
figure(1611);clf
    compare(sysID.data.CL.valid.direct,G_CL_dir)
figure(1612);clf
    resid(sysID.data.CL.valid.direct,G_CL_dir)

%% for indirect method prep

% Make the data for identification from r to y (input to the loop to output of the loop)
%   This results in the process sensitivity of the system
sysID.data.CL.train.indirect1 = iddata(sysID.data.CL.train.y,sysID.data.CL.train.r,Ts);
sysID.data.CL.valid.indirect1 = iddata(sysID.data.CL.valid.y,sysID.data.CL.valid.r,Ts);


% Perform a parametric identification with an OE model structure, since for
%   the indirect methods the noise model do not matter anymore therefore they
%   can be neglected in the model structure as well

% est_yr = oe(dataCL_direct_yr, [nb nf 0]);
est_yr = oe(sysID.data.CL.train.indirect1, [3 3 0]);
G_yr = tf(est_yr.B,est_yr.F,est_yr.Ts,'variable','q^-1');

% Verification using the residial test, only look at the cross correlation
figure(1613);
    resid(sysID.data.CL.valid.indirect1, est_yr)

%% for indirect method prep

% Make the data for identification from r to u (input to the loop to input of the plant)
%   This results in the input sensitivity of the system
sysID.data.CL.train.indirect2 = iddata(sysID.data.CL.train.u,sysID.data.CL.train.r,Ts);
sysID.data.CL.valid.indirect2 = iddata(sysID.data.CL.valid.u,sysID.data.CL.valid.r,Ts);


% Perform a parametric identification with an OE model structure, since for
%   the indirect methods the noise model do not matter anymore therefore they
%   can be neglected in the model structure as well

% est_ur = oe(dataCL_direct_ur, [nb nf 0]);
est_ur = oe(sysID.data.CL.train.indirect2, [4 4 0]);

G_ur = tf(est_ur.B,est_ur.F,est_ur.Ts,'variable','q^-1');

% Verification using the residial test, only look at the cross correlation
figure(1614);
    resid(sysID.data.CL.valid.indirect2, est_ur)

%% Indirect approach 1 - co-prime
clc
% Since G_yr is the process sensitivity and G_ur is the input sensitivity;
%   the division between those two results in the plant model, this is the
%   essence of the co-prime method
G_CL_cp = G_yr/G_ur;

% Verify with previously obtainted model
figure(1615);
    bodemag(G_CL_cp,G(1,1));
    legend('Co-prime, closed loop','Non-param, open loop')
figure(1616);clf
    compare(sysID.data.CL.valid.direct,G_CL_cp)
figure(2616);clf
    resid(sysID.data.CL.valid.direct,G_CL_cp)

%% Indirect approach 2 - two-stage

% Simulate the plant input (instead of using the measurement) using the 
%   input sensitivity as defined before
u_r = lsim(G_ur,sysID.data.CL.train.r,round(sysID.data.CL.train.tVec,3));

% Use the simulated plant input data for the identification in combination
%   for the measured output, known as the two-stage method
dataCL_direct_2s = iddata(sysID.data.CL.train.y,u_r,Ts);

% Only a OE model structure should be enough since the noise model is not
%   needed for the indirect identification methods
est_2s = oe(dataCL_direct_2s, [2 4 0]);

G_CL_2s = tf(est_2s.B,est_2s.F,est_2s.Ts);

% Verify with previously obtainted model
figure(617);clf;hold on;
    bode(G_CL_2s,G(1,1),opt)
    xlim([bodeRange(1) bodeRange(end)])
figure(1617);clf
    compare(sysID.data.CL.valid.direct,G_CL_2s)
figure(2617);clf
    resid(sysID.data.CL.valid.direct,G_CL_2s)


%%
figure(1001);clf
    bodemag(G_CL_dir);hold on;grid minor
    bodemag(G_CL_cp)
    bodemag(G_CL_2s)
    bodemag(G(1,1),'g--')
        xlim([1e-6 1e0])
        legend('CL, direct','CL, co-prime','CL, 2-stage','Model','location','south')

figure(1002);clf
    compare(sysID.data.CL.train.direct,G_CL_dir,G_CL_cp,G_CL_2s,G(1,1));hold on;grid minor
        legend('Org. data','CL, direct','CL, co-prime','CL, 2-stage','Model','location','south')

%% CHECK SLIDES
% lecture 11, slide 26























































