
%% Import data and format for later use
% dir1 = 'p__CB_none__ST_clms__SM_2w__SA_1w__DT_251009__MD_8h__WT_no__DS_10';
dir1 = 'p__CB_none__ST_clms__SM_2w__SA_1w__DT_251010__MD_8h__WT_no__DS_10';
dir2 = 'C:\Users\maxja\Documents\(4)School\Master\Q9_Internship\MaxJans_InternshipCS_2025\SystemIdentification\_Data\9Okt25\usedController.mat';
[sysID.CLdata, sysID.CLdataFilt, C] = makeCLmeasData(dir1,dir2);



%% Do a non-parametric identification
[sysID.nonPar.CL.frf.direct.raw,~]  = makeClosedLoopFRF(sysID.CLdata.train.e,sysID.CLdata.train.r,Ck.tot,fs);
[sysID.nonPar.CL.frf.classic.raw,~] = makeOpenLoopFRF_sysIdent(sysID.CLdata.train.y,sysID.CLdata.train.u,fs);
[G_coprime_ru_raw,~]                = makeOpenLoopFRF_sysIdent(sysID.CLdata.train.u,sysID.CLdata.train.r,fs);
[G_coprime_ry_raw,~]                = makeOpenLoopFRF_sysIdent(sysID.CLdata.train.y,sysID.CLdata.train.r,fs);
sysID.nonPar.CL.frf.coprime.raw     = G_coprime_ry_raw/G_coprime_ru_raw;

[sysID.nonPar.CL.frf.direct.filt,~]  = makeClosedLoopFRF(sysID.CLdataFilt.train.e,sysID.CLdataFilt.train.r,Ck.tot,fs);
[sysID.nonPar.CL.frf.classic.filt,~] = makeOpenLoopFRF_sysIdent(sysID.CLdataFilt.train.y,sysID.CLdataFilt.train.u,fs);
[G_coprime_ru,~]                     = makeOpenLoopFRF_sysIdent(sysID.CLdataFilt.train.u,sysID.CLdataFilt.train.r,fs);
[G_coprime_ry,~]                     = makeOpenLoopFRF_sysIdent(sysID.CLdataFilt.train.y,sysID.CLdataFilt.train.r,fs);
sysID.nonPar.CL.frf.coprime.filt     = G_coprime_ry/G_coprime_ru;

if plotAll 
    CL_sim_validate_nonParam
end

%% Do a parametric identification, direct method
sysID.CLdata.train.direct     = iddata(sysID.CLdata.train.y,sysID.CLdata.train.u,Ts);
sysID.CLdata.valid.direct     = iddata(sysID.CLdata.valid.y,sysID.CLdata.valid.u,Ts);
sysID.CLdata.trans.direct     = iddata(sysID.CLdata.trans.y,sysID.CLdata.trans.u,Ts);

sysID.CLdataFilt.train.direct = iddata(sysID.CLdataFilt.train.y,sysID.CLdataFilt.train.u,Ts);
sysID.CLdataFilt.valid.direct = iddata(sysID.CLdataFilt.valid.y,sysID.CLdataFilt.valid.u,Ts);
sysID.CLdataFilt.trans.direct = iddata(sysID.CLdataFilt.trans.y,sysID.CLdataFilt.trans.u,Ts);

% structure
nb = 1;
nc = 2;
nd = 1;
nf = 3;
nk = 0;
est_raw = bj(sysID.CLdata.train.direct, [nb nc nd nf nk]);
est = bj(sysID.CLdataFilt.train.direct, [nb nc nd nf nk]);

% Make the transfer function with the identified parameters
sysID.par.CL.time.direct.raw  = tf(est_raw.B,est_raw.F,est_raw.Ts,'variable','q^-1');
sysID.par.CL.time.direct.filt = tf(est.B,est.F,est.Ts,'variable','q^-1');

H_CL_dir = tf(est.C,est.D,est.Ts,'variable','q^-1');

if plotAll 
    CL_sim_validate_direct_param
end

% figure(1);clf; resid(sysID.CLdata.valid.direct,    sysID.par.CL.time.direct.raw)
% figure(2);clf; resid(sysID.CLdataFilt.valid.direct,sysID.par.CL.time.direct.filt)

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
u_r_raw = lsim(G_ur_raw,sysID.CLdata.train.r,round(sysID.CLdata.train.tVec,3));
u_r     = lsim(G_ur,sysID.CLdataFilt.train.r,round(sysID.CLdataFilt.train.tVec,3));

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





