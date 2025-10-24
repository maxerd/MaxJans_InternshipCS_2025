function [data, dataFilt, Ck, varargout] = makeCLmeasData(dataDir_CL,contDir)

%% Define the path of the controller
Ck = load(contDir);

%% Some user inputs
makeValSet  = true;  % [true/false] Is a seperate validation set needed?
removeTrans = true; % [true/false] Does the transient need to be removed

%% Define the measurement paths
% Load in the measurement data
disp('Loading in measurement data...')
    dat  = load(dataDir_CL);
disp('Loading done!')

%% Define some simulation variables

T0 = 23; % [degC] Inital Temperature
Tamb = 23; % [degC] Ambient Temperature

%% Define some measurement variables for easier use
Ts = round(dat.Ts,0);  % [s]     Sampling time
fs = 1/Ts;     % [Hz]    Sampling frequency
N = size(dat.tVec,2); % [-] Measurement samples

tMeas = N*Ts/3600; % [hours] Measurement time

tVec = dat.tVec; % [s] Time vector

dist = dat.dist+55;

%% Simulate the systems

% Ambient temperature vector, constant ambient temperature for simulations
    ambVec = dat.tempAmb;

% Run the closed loop simulation, again with the previously defined noise signal added
    y_CL   = dat.tempTM;

%% Define the in and output vector sizes

% Define the range of the data that is used for identification and validation
if removeTrans
    idx1_CL = round(1.25*3600*fs,0);
    % idx1_CL = round(0.5*3600*fs,0);
else
    idx1_CL = 1;
end
if makeValSet
    idx2_CL = round(N-(N-2.5*3600*fs)*0.2,0);
else
    idx2_CL = N;
end

idxRange_CL_trns = 1:idx1_CL+1;
idxRange_CL      = idx1_CL:idx2_CL;
idxRange_CL_val  = (idx2_CL):N;

%% Define the in and output vectors, for easy and consistent use
data.trans.y = (y_CL(idxRange_CL_trns)-ambVec(idxRange_CL_trns))';
data.trans.u = (dat.Watt(idxRange_CL_trns))';
data.trans.d = (dist(idxRange_CL_trns))';
data.trans.e = (dat.error(idxRange_CL_trns))';
data.trans.r = (data.trans.e+data.trans.y);

data.train.y = (y_CL(idxRange_CL)-ambVec(idxRange_CL))';
data.train.u = (dat.Watt(idxRange_CL))';
data.train.d = (dist(idxRange_CL))';
data.train.e = (dat.error(idxRange_CL))';
data.train.r = (data.train.e+data.train.y);

data.valid.y = (y_CL(idxRange_CL_val)-ambVec(idxRange_CL_val))';
data.valid.u = (dat.Watt(idxRange_CL_val))';
data.valid.d = (dist(idxRange_CL_val))';
data.valid.e = (dat.error(idxRange_CL_val))';
data.valid.r = (data.valid.e+data.valid.y);

data.full.y = (y_CL-ambVec)';
data.full.u = (dat.Watt)';
data.full.d = (dist)';
data.full.e = (dat.error)';
data.full.r = (data.full.e+data.full.y);

% Define the time vectors
data.trans.tVec = tVec(idxRange_CL_trns);
data.train.tVec = tVec(idxRange_CL);
data.valid.tVec = tVec(idxRange_CL_val);
data.full.tVec  = tVec;

%% Filter the data vectors to remove some of the effect of noise
% High order lowpass filter, to filter away all the noise above the 
% multisine exitation frequency
    LPfilt = orderLP(0.01,2,0,0,0.7)^2;

% Actually filter the data
    dataFilt.trans.y = lsim(LPfilt,data.trans.y-data.trans.y(1),data.trans.tVec)+data.trans.y(1);
    dataFilt.trans.u = lsim(LPfilt,data.trans.u-data.trans.u(1),data.trans.tVec)+data.trans.u(1);
    dataFilt.trans.d = lsim(LPfilt,data.trans.d-data.trans.d(1),data.trans.tVec)+data.trans.d(1);
    dataFilt.trans.e = lsim(LPfilt,data.trans.e-data.trans.e(1),data.trans.tVec)+data.trans.e(1);
    dataFilt.trans.r = lsim(LPfilt,data.trans.r-data.trans.r(1),data.trans.tVec)+data.trans.r(1);

    dataFilt.train.y = lsim(LPfilt,data.train.y-data.train.y(1),data.train.tVec)+data.train.y(1);
    dataFilt.train.u = lsim(LPfilt,data.train.u-data.train.u(1),data.train.tVec)+data.train.u(1);
    dataFilt.train.d = lsim(LPfilt,data.train.d-data.train.d(1),data.train.tVec)+data.train.d(1);
    dataFilt.train.e = lsim(LPfilt,data.train.e-data.train.e(1),data.train.tVec)+data.train.e(1);
    dataFilt.train.r = lsim(LPfilt,data.train.r-data.train.r(1),data.train.tVec)+data.train.r(1);

    dataFilt.valid.y = lsim(LPfilt,data.valid.y-data.valid.y(1),data.valid.tVec)+data.valid.y(1);
    dataFilt.valid.u = lsim(LPfilt,data.valid.u-data.valid.u(1),data.valid.tVec)+data.valid.u(1);
    dataFilt.valid.d = lsim(LPfilt,data.valid.d-data.valid.d(1),data.valid.tVec)+data.valid.d(1);
    dataFilt.valid.e = lsim(LPfilt,data.valid.e-data.valid.e(1),data.valid.tVec)+data.valid.e(1);
    dataFilt.valid.r = lsim(LPfilt,data.valid.r-data.valid.r(1),data.valid.tVec)+data.valid.r(1);

    dataFilt.full.y = lsim(LPfilt,data.full.y-data.full.y(1),data.full.tVec)+data.full.y(1);
    dataFilt.full.u = lsim(LPfilt,data.full.u-data.full.u(1),data.full.tVec)+data.full.u(1);
    dataFilt.full.d = lsim(LPfilt,data.full.d-data.full.d(1),data.full.tVec)+data.full.d(1);
    dataFilt.full.e = lsim(LPfilt,data.full.e-data.full.e(1),data.full.tVec)+data.full.e(1);
    dataFilt.full.r = lsim(LPfilt,data.full.r-data.full.r(1),data.full.tVec)+data.full.r(1);
    
% Define the time vectors
dataFilt.trans.tVec = data.trans.tVec;
dataFilt.train.tVec = data.train.tVec;
dataFilt.valid.tVec = data.valid.tVec;
dataFilt.full.tVec  = data.full.tVec ;


































end