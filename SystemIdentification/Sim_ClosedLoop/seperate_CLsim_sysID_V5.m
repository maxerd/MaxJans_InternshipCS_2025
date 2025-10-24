%——————————————————————————————————————————————————————————————————————————————————————————
%  ██████╗██╗      ██████╗  ██████╗███████╗██████╗      ██╗      ██████╗  ██████╗ ██████╗ 
% ██╔════╝██║     ██╔═══██╗██╔════╝██╔════╝██╔══██╗     ██║     ██╔═══██╗██╔═══██╗██╔══██╗
% ██║     ██║     ██║   ██║╚█████╗ ███████╗██║  ██║     ██║     ██║   ██║██║   ██║██████╔╝
% ██║     ██║     ██║   ██║ ╚═══██╗██╔════╝██║  ██║     ██║     ██║   ██║██║   ██║██╔═══╝ 
% ╚██████╗███████╗╚██████╔╝██████╔╝███████╗██████╔╝     ███████╗╚██████╔╝╚██████╔╝██║     
%  ╚═════╝╚══════╝ ╚═════╝ ╚═════╝ ╚══════╝╚═════╝      ╚══════╝ ╚═════╝  ╚═════╝ ╚═╝     
%
%  ██████╗██╗███╗   ███╗██╗   ██╗██╗      █████╗ ████████╗██╗ ██████╗ ███╗  ██╗
% ██╔════╝██║████╗ ████║██║   ██║██║     ██╔══██╗╚══██╔══╝██║██╔═══██╗████╗ ██║
% ╚█████╗ ██║██╔████╔██║██║   ██║██║     ██║  ██║   ██║   ██║██║   ██║██╔██╗██║
%  ╚═══██╗██║██║╚██╔╝██║██║   ██║██║     ███████║   ██║   ██║██║   ██║██║╚████║
% ██████╔╝██║██║ ╚═╝ ██║╚██████╔╝███████╗██╔══██║   ██║   ██║╚██████╔╝██║ ╚███║
% ╚═════╝ ╚═╝╚═╝     ╚═╝ ╚═════╝ ╚══════╝╚═╝  ╚═╝   ╚═╝   ╚═╝ ╚═════╝ ╚═╝  ╚══╝
%——————————————————————————————————————————————————————————————————————————————————————————

%% Run simulation and format data
% Always uses the pre-defined disturbance and noise signal, since these
% might be needed to stay constant over several simulations

% Define the bandwidth of the controller that is going to be generated
controllerBW = 0.001;

% Generate the data itself
if newData
    [sysID.CLdata, sysID.CLdataFilt, C] = makeCLdata(G,tMeas,controllerBW,dist,v,Ts);
end

%% Do a non-parametric identification on non-transient-compensated data
% Perform non-parametric identifications on the raw non-compensated data
    % Direct
    [sysID.nonPar.CL.frf.trd.direct.raw ,~] = makeOpenLoopFRF_sysIdent(sysID.CLdata.train.y,sysID.CLdata.train.u,fs);
    
    % Indirect, classical
    [sysID.nonPar.CL.frf.trd.classic.raw,~] = makeClosedLoopFRF(sysID.CLdata.train.e,sysID.CLdata.train.r,Ck.tot,fs); 

    % Indirect, co-prime
    [G_coprime_ru_raw               ,~] = makeOpenLoopFRF_sysIdent(sysID.CLdata.train.u,sysID.CLdata.train.r,fs); 
    [G_coprime_ry_raw               ,~] = makeOpenLoopFRF_sysIdent(sysID.CLdata.train.y,sysID.CLdata.train.r,fs);
    sysID.nonPar.CL.frf.trd.coprime.raw = G_coprime_ry_raw/G_coprime_ru_raw;

    % LPM versions
    polyOrder = 6; % The order of the polynomial that is fitted to the data
    locality  = 8; % Amount of points (pos & neg) to consider around apprx freq

    % Direct
    [sysID.nonPar.CL.frf.lpm.direct.raw,~] = sysID_LPM(sysID.CLdata.train.y,sysID.CLdata.train.u',fs,locality,polyOrder);
    
    % Indirect, classical
    [G_classic_du_raw,~]                 = sysID_LPM(sysID.CLdata.train.e,sysID.CLdata.train.r',fs,locality,polyOrder);
    sysID.nonPar.CL.frf.lpm.classic.raw   = (inv(G_classic_du_raw)-1)*inv(Ck.tot);

    % Indirect, co-prime
    [G_coprime_ru_raw,~]                 = sysID_LPM(sysID.CLdata.train.u,sysID.CLdata.train.r',fs,locality,polyOrder);
    [G_coprime_ry_raw,~]                 = sysID_LPM(sysID.CLdata.train.y,sysID.CLdata.train.r',fs,locality,polyOrder);
    sysID.nonPar.CL.frf.lpm.coprime.raw  = G_coprime_ry_raw/G_coprime_ru_raw;

% Perform non-parametric identifications on the filtered non-compensated data
    % Direct
    [sysID.nonPar.CL.frf.trd.direct.filt,~] = makeOpenLoopFRF_sysIdent(sysID.CLdataFilt.train.y,sysID.CLdataFilt.train.u,fs);
    
    % Indirect, classical
    [sysID.nonPar.CL.frf.trd.classic.filt,~]  = makeClosedLoopFRF(sysID.CLdataFilt.train.e,sysID.CLdataFilt.train.r,Ck.tot,fs);

    % Indirect, co-prime
    [G_coprime_ru,~]                     = makeOpenLoopFRF_sysIdent(sysID.CLdataFilt.train.u,sysID.CLdataFilt.train.r,fs);
    [G_coprime_ry,~]                     = makeOpenLoopFRF_sysIdent(sysID.CLdataFilt.train.y,sysID.CLdataFilt.train.r,fs);
    sysID.nonPar.CL.frf.trd.coprime.filt     = G_coprime_ry/G_coprime_ru;

    % Direct
    [sysID.nonPar.CL.frf.lpm.direct.filt,~] = makeOpenLoopFRF_sysIdent(sysID.CLdataFilt.train.y,sysID.CLdataFilt.train.u,fs);

    % Indirect, classical
    [G_classic_du,~]                 = sysID_LPM(sysID.CLdataFilt.train.e,sysID.CLdataFilt.train.r',fs,locality,polyOrder);
    sysID.nonPar.CL.frf.lpm.classic.filt   = (Ck.tot^(-1))*((G_classic_du^(-1))-1);

    % Indirect, co-prime
    [G_coprime_ru,~]                     = sysID_LPM(sysID.CLdataFilt.train.u,sysID.CLdataFilt.train.r',fs,locality,polyOrder);
    [G_coprime_ry,~]                     = sysID_LPM(sysID.CLdataFilt.train.y,sysID.CLdataFilt.train.r',fs,locality,polyOrder);
    sysID.nonPar.CL.frf.lpm.coprime.filt = G_coprime_ry/G_coprime_ru;

%% Compensate for first order behaviour
if compTrans
    % Make the complementary, control and input sensitivity functions (respectively)
    CLsys = feedback(sysID.par.OL.time.firstApprox.raw*Ck.tot,1);
    CL_CS = feedback(Ck.tot,sysID.par.OL.time.firstApprox.raw);
    CL_S  = feedback(1,sysID.par.OL.time.firstApprox.raw*Ck.tot);

    % Simulate the complementary, control and input sensitivity
    yComp = lsim(CLsys,(mean(sysID.CLdata.train.r)).*ones(size(sysID.CLdata.train.tVec)),sysID.CLdata.train.tVec);
    uComp = lsim(CL_CS,(mean(sysID.CLdata.train.r)).*ones(size(sysID.CLdata.train.tVec)),sysID.CLdata.train.tVec);
    eComp = lsim(CL_S,(mean(sysID.CLdata.train.r)).*ones(size(sysID.CLdata.train.tVec)),sysID.CLdata.train.tVec);
    
    % Calculate the steady state value of the signals that are going to be compensated
    yDC = mean(sysID.CLdata.train.y(end-1000:end));
    uDC = mean(sysID.CLdata.train.u(end-1000:end));
    eDC = mean(sysID.CLdata.train.e(end-1000:end));
    
    % Substract the first order behavior from the actual (simulation) data
    sysID.CLdata.train.y     = sysID.CLdata.train.y+yDC-yComp;
    sysID.CLdataFilt.train.y = sysID.CLdataFilt.train.y+yDC-yComp;
    sysID.CLdata.train.u     = sysID.CLdata.train.u+uDC-uComp;
    sysID.CLdataFilt.train.u = sysID.CLdataFilt.train.u+uDC-uComp;
    sysID.CLdata.train.e     = sysID.CLdata.train.e+eDC-yComp;
    sysID.CLdataFilt.train.e = sysID.CLdataFilt.train.e+eDC-yComp;
end

%% Do a non-parametric identification on transient-compensated data
% Perform non-parametric identifications on the raw non-compensated data
    % Direct
    [sysID.nonPar.CL.frfComp.trd.direct.raw ,~] = makeOpenLoopFRF_sysIdent(sysID.CLdata.train.y,sysID.CLdata.train.u,fs);
    
    % Indirect, classical
    [sysID.nonPar.CL.frfComp.trd.classic.raw,~] = makeClosedLoopFRF(sysID.CLdata.train.e,sysID.CLdata.train.r,Ck.tot,fs); 

    % Indirect, co-prime
    [G_coprime_ru_raw               ,~] = makeOpenLoopFRF_sysIdent(sysID.CLdata.train.u,sysID.CLdata.train.r,fs); 
    [G_coprime_ry_raw               ,~] = makeOpenLoopFRF_sysIdent(sysID.CLdata.train.y,sysID.CLdata.train.r,fs);
    sysID.nonPar.CL.frfComp.trd.coprime.raw = G_coprime_ry_raw/G_coprime_ru_raw;

    % LPM versions
    polyOrder = 6; % The order of the polynomial that is fitted to the data
    locality  = 8; % Amount of points (pos & neg) to consider around apprx freq

    % Direct
    [sysID.nonPar.CL.frfComp.lpm.direct.raw,~] = sysID_LPM(sysID.CLdata.train.y,sysID.CLdata.train.u',fs,locality,polyOrder);
    
    % Indirect, classical
    [G_classic_du_raw,~]                 = sysID_LPM(sysID.CLdata.train.e,sysID.CLdata.train.r',fs,locality,polyOrder);
    sysID.nonPar.CL.frfComp.lpm.classic.raw   = (inv(G_classic_du_raw)-1)*inv(Ck.tot);

    % Indirect, co-prime
    [G_coprime_ru_raw,~]                 = sysID_LPM(sysID.CLdata.train.u,sysID.CLdata.train.r',fs,locality,polyOrder);
    [G_coprime_ry_raw,~]                 = sysID_LPM(sysID.CLdata.train.y,sysID.CLdata.train.r',fs,locality,polyOrder);
    sysID.nonPar.CL.frfComp.lpm.coprime.raw  = G_coprime_ry_raw/G_coprime_ru_raw;

% Perform non-parametric identifications on the filtered non-compensated data
    % Direct
    [sysID.nonPar.CL.frfComp.trd.direct.filt,~] = makeOpenLoopFRF_sysIdent(sysID.CLdataFilt.train.y,sysID.CLdataFilt.train.u,fs);
    
    % Indirect, classical
    [sysID.nonPar.CL.frfComp.trd.classic.filt,~]  = makeClosedLoopFRF(sysID.CLdataFilt.train.e,sysID.CLdataFilt.train.r,Ck.tot,fs);

    % Indirect, co-prime
    [G_coprime_ru,~]                     = makeOpenLoopFRF_sysIdent(sysID.CLdataFilt.train.u,sysID.CLdataFilt.train.r,fs);
    [G_coprime_ry,~]                     = makeOpenLoopFRF_sysIdent(sysID.CLdataFilt.train.y,sysID.CLdataFilt.train.r,fs);
    sysID.nonPar.CL.frfComp.trd.coprime.filt     = G_coprime_ry/G_coprime_ru;

    % Direct
    [sysID.nonPar.CL.frfComp.lpm.direct.filt,~] = makeOpenLoopFRF_sysIdent(sysID.CLdataFilt.train.y,sysID.CLdataFilt.train.u,fs);

    % Indirect, classical
    [G_classic_du,~]                 = sysID_LPM(sysID.CLdataFilt.train.e,sysID.CLdataFilt.train.r',fs,locality,polyOrder);
    sysID.nonPar.CL.frfComp.lpm.classic.filt   = (Ck.tot^(-1))*((G_classic_du^(-1))-1);

    % Indirect, co-prime
    [G_coprime_ru,~]                     = sysID_LPM(sysID.CLdataFilt.train.u,sysID.CLdataFilt.train.r',fs,locality,polyOrder);
    [G_coprime_ry,~]                     = sysID_LPM(sysID.CLdataFilt.train.y,sysID.CLdataFilt.train.r',fs,locality,polyOrder);
    sysID.nonPar.CL.frfComp.lpm.coprime.filt = G_coprime_ry/G_coprime_ru;

%% Do a parametric identification, direct method
% Define the iddata datasets for the direct identification method
sysID.CLdata.train.direct     = iddata(sysID.CLdata.train.y,sysID.CLdata.train.u,Ts);
sysID.CLdata.valid.direct     = iddata(sysID.CLdata.valid.y,sysID.CLdata.valid.u,Ts);
sysID.CLdata.trans.direct     = iddata(sysID.CLdata.trans.y,sysID.CLdata.trans.u,Ts);

sysID.CLdataFilt.train.direct = iddata(sysID.CLdataFilt.train.y,sysID.CLdataFilt.train.u,Ts);
sysID.CLdataFilt.valid.direct = iddata(sysID.CLdataFilt.valid.y,sysID.CLdataFilt.valid.u,Ts);
sysID.CLdataFilt.trans.direct = iddata(sysID.CLdataFilt.trans.y,sysID.CLdataFilt.trans.u,Ts);

% Box-jenkins structure with 
nb = 3;
nc = 2;
nd = 1;
nf = 3;
nk = 0;
est_raw = bj(sysID.CLdata.train.direct, [nb nc nd nf nk]);
est = bj(sysID.CLdataFilt.train.direct, [nb nc nd nf nk]);

% Make the transfer function with the identified parameters
sysID.par.CL.time.direct.raw  = tf(est_raw.B,est_raw.F,est_raw.Ts,'variable','q^-1');
sysID.par.CL.time.direct.filt = tf(est.B,est.F,est.Ts,'variable','q^-1');

sysID.par.CL.time.direct.raw

H_CL_dir = tf(est.C,est.D,est.Ts,'variable','q^-1');

if plotAll 
    CL_sim_validate_direct_param
end

figure(1);clf; compare(sysID.CLdata.valid.direct,    sysID.par.CL.time.direct.raw,sysID.par.CL.time.direct.filt,1)
figure(2);clf; compare(sysID.CLdataFilt.valid.direct,sysID.par.CL.time.direct.raw,sysID.par.CL.time.direct.filt,1)

%% for indirect method prep

% Make the data for identification from r to y (input to the loop to output of the loop)
%   This results in the process sensitivity of the system
sysID.CLdata.train.indirect1     = iddata(sysID.CLdata.train.y,sysID.CLdata.train.r,Ts);
sysID.CLdata.valid.indirect1     = iddata(sysID.CLdata.valid.y,sysID.CLdata.valid.r,Ts);

sysID.CLdataFilt.train.indirect1 = iddata(sysID.CLdataFilt.train.y,sysID.CLdataFilt.train.r,Ts);
sysID.CLdataFilt.valid.indirect1 = iddata(sysID.CLdataFilt.valid.y,sysID.CLdataFilt.valid.r,Ts);


% Perform a parametric identification with an OE model structure, since for
%   the indirect methods the noise model do not matter anymore therefore they
%   can be neglected in the model structure as well

est_yr_raw = oe(sysID.CLdata.train.indirect1, [3 3 0]);
G_yr_raw   = tf(est_yr_raw.B,est_yr_raw.F,est_yr_raw.Ts,'variable','q^-1');

est_yr = oe(sysID.CLdata.train.indirect1, [3 3 0]);
G_yr   = tf(est_yr.B,est_yr.F,est_yr.Ts,'variable','q^-1');

% figure(1);clf; resid(sysID.CLdata.valid.indirect1,G_yr_raw)
% figure(2);clf; resid(sysID.CLdataFilt.valid.indirect1,G_yr)

%% for indirect method prep

% Make the data for identification from r to u (input to the loop to input of the plant)
%   This results in the input sensitivity of the system
sysID.CLdata.train.indirect2     = iddata(sysID.CLdata.train.u,sysID.CLdata.train.r,Ts);
sysID.CLdata.valid.indirect2     = iddata(sysID.CLdata.valid.u,sysID.CLdata.valid.r,Ts);

sysID.CLdataFilt.train.indirect2 = iddata(sysID.CLdataFilt.train.u,sysID.CLdataFilt.train.r,Ts);
sysID.CLdataFilt.valid.indirect2 = iddata(sysID.CLdataFilt.valid.u,sysID.CLdataFilt.valid.r,Ts);


% Perform a parametric identification with an OE model structure, since for
%   the indirect methods the noise model do not matter anymore therefore they
%   can be neglected in the model structure as well

est_ur_raw = oe(sysID.CLdata.train.indirect2, [2 3 0]);
G_ur_raw   = tf(est_ur_raw.B,est_ur_raw.F,est_ur_raw.Ts,'variable','q^-1');

est_ur = oe(sysID.CLdataFilt.train.indirect2, [3 3 0]);
G_ur   = tf(est_ur.B,est_ur.F,est_ur.Ts,'variable','q^-1');

% figure(1);clf; resid(sysID.CLdata.valid.indirect2,G_ur_raw)
% figure(2);clf; resid(sysID.CLdataFilt.valid.indirect2,G_ur)

%% Indirect approach 3 - classical
sysID.CLdata.train.classic     = iddata(sysID.CLdata.train.u,sysID.CLdata.train.d,Ts);
sysID.CLdata.valid.classic     = iddata(sysID.CLdata.valid.u,sysID.CLdata.valid.d,Ts);

sysID.CLdataFilt.train.classic = iddata(sysID.CLdataFilt.train.u,sysID.CLdataFilt.train.d,Ts);
sysID.CLdataFilt.valid.classic = iddata(sysID.CLdataFilt.valid.u,sysID.CLdataFilt.valid.d,Ts);

est_du_raw = oe(sysID.CLdata.train.classic, [5 6 0]);
G_du_raw   = tf(est_du_raw.B,est_du_raw.F,est_du_raw.Ts,'variable','q^-1');

est_du = oe(sysID.CLdata.train.classic, [5 6 0]);
G_du   = tf(est_du.B,est_du.F,est_du.Ts,'variable','q^-1');


C_classic = c2d(Ck.tot,Ts);

sysID.par.CL.time.classic.raw = (1-inv(G_du_raw))*inv(C_classic);
sysID.par.CL.time.classic.filt  = (1-inv(G_du))*inv(C_classic);

% figure(1);clf; resid(sysID.CLdata.valid.classic,G_du_raw)
% figure(2);clf; resid(sysID.CLdataFilt.valid.classic,G_du)
% figure(1);clf; resid(sysID.CLdata.valid.direct,    sysID.par.CL.time.classic.raw)
% figure(2);clf; resid(sysID.CLdataFilt.valid.direct,sysID.par.CL.time.classic.filt)

%% Indirect approach 1 - co-prime
% Since G_yr is the process sensitivity and G_ur is the input sensitivity;
%   the division between those two results in the plant model, this is the
%   essence of the co-prime method
sysID.par.CL.time.coprime.raw = G_yr_raw/G_ur_raw;
sysID.par.CL.time.coprime.filt = G_yr/G_ur;

if plotAll 
    CL_sim_validate_coprime_param
end

% figure(1);clf; resid(sysID.CLdata.valid.direct,    sysID.par.CL.time.coprime.raw)
% figure(2);clf; resid(sysID.CLdataFilt.valid.direct,sysID.par.CL.time.coprime.filt)

%% Indirect approach 2 - two-stage

% Simulate the plant input (instead of using the measurement) using the 
%   input sensitivity as defined before
u_r_raw = lsim(G_ur_raw,sysID.CLdata.train.r,round(sysID.CLdata.train.tVec,0));
u_r     = lsim(G_ur,sysID.CLdataFilt.train.r,round(sysID.CLdataFilt.train.tVec,0));

% Use the simulated plant input data for the identification in combination
%   for the measured output, known as the two-stage method
dataCL_direct_2s_raw = iddata(sysID.CLdata.train.y,u_r_raw,Ts);
dataCL_direct_2s     = iddata(sysID.CLdataFilt.train.y,u_r,Ts);

% Only a OE model structure should be enough since the noise model is not
%   needed for the indirect identification methods
est_2s_raw = oe(dataCL_direct_2s_raw, [2 4 0]);
est_2s     = oe(dataCL_direct_2s, [2 4 0]);

sysID.par.CL.time.twostage.raw = tf(est_2s_raw.B,est_2s_raw.F,est_2s_raw.Ts);
sysID.par.CL.time.twostage.filt = tf(est_2s.B,est_2s.F,est_2s.Ts);

if plotAll 
    CL_sim_validate_twoStage_param
end


% figure(1);clf; resid(dataCL_direct_2s_raw,est_2s_raw)
% figure(2);clf; resid(dataCL_direct_2s,    est_2s)
% figure(1);clf; resid(sysID.CLdata.valid.direct,    sysID.par.CL.time.twostage.raw)
% figure(2);clf; resid(sysID.CLdataFilt.valid.direct,sysID.par.CL.time.twostage.filt)





