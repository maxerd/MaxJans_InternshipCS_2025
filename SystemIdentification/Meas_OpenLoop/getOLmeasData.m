function [data, dataFilt] = getOLmeasData(dataDir_OL)

%% Some user inputs
makeValSet  = true;  % [true/false] Is a seperate validation set needed?
removeTrans = false; % [true/false] Does the transient need to be removed


%% Define the measurement paths
% dataDir_OL = 'C:\Users\maxja\Documents\(4)School\Master\Q9_Internship\matlabFiles\measurements\processedData\p__CB_none__ST_ms__SM_24w__SA_12w__DT_250829__MD_8h__WT_no__DS_1000.mat'; % Open loop data used for identification

% Load in the measurement data
disp('Loading in measurement data...')
    dat   = load(dataDir_OL);
disp('Loading done!')

%% Define some measurement variables for easier use
Ts = round(dat.Ts,0);  % [s]     Sampling time
fs = 1/Ts;     % [Hz]    Sampling frequency
N = size(dat.tVec,2); % [-] Measurement samples

tMeas = N*Ts/3600; % [hours] Measurement time

tVec = dat.tVec; % [s] Time vector

%% Define the input signal
% dist     = dat.Watt./5;
dist     = dat.Watt;

%% Define some data vectors
% To correspond to the simulation file as well

% Ambient temperature vector
    ambVec = dat.tempAmb;

% Define the system output, for easier use later
    y_OL   = dat.tempTM;

%% Define the range of the data that is used for identification and validation
if removeTrans
    idx1_OL = 8*3600*fs;
else
    idx1_OL = 1;
end
if makeValSet
    idx2_OL = round(N-(N-6*3600*fs)*0.2);
else
    idx2_OL = N;
end

idxRange_OL_trns = 1:idx1_OL+1;     % Transient data range
idxRange_OL      = idx1_OL:idx2_OL; % Training data range
idxRange_OL_val  = (idx2_OL-1):N;   % Validation data range

%% Define the in and output vectors, for easy and consistent use
data.trans.y = (y_OL(idxRange_OL_trns)-ambVec(idxRange_OL_trns))';
data.trans.u  = (dist(idxRange_OL_trns))';

data.train.y = (y_OL(idxRange_OL)-ambVec(idxRange_OL))';
data.train.u  = (dist(idxRange_OL))';

data.valid.y = (y_OL(idxRange_OL_val)-ambVec(idxRange_OL_val))';
data.valid.u  = (dist(idxRange_OL_val))';

data.full.y = (y_OL-ambVec)';
data.full.u  = (dist)';

% Define the time vectors
data.trans.tVec = tVec(idxRange_OL_trns);
data.train.tVec = tVec(idxRange_OL);
data.valid.tVec   = tVec(idxRange_OL_val);
data.full.tVec  = tVec;

%% Filter the data vectors to remove some of the effect of noise
% High order lowpass filter, to filter away all the noise above the 
% multisine exitation frequency
    LPfilt = orderLP(0.01,2,0,0,0.7)^2;

% Actually filter the data
    % Open loop data
    dataFilt.trans.y = lsim(LPfilt,data.trans.y-data.trans.y(1),data.trans.tVec)+data.trans.y(1);
    dataFilt.trans.u  = lsim(LPfilt,data.trans.u-data.trans.u(1)  ,data.trans.tVec)+data.trans.u(1);

    dataFilt.train.y = lsim(LPfilt,data.train.y-data.train.y(1),data.train.tVec)+data.train.y(1);
    dataFilt.train.u  = lsim(LPfilt,data.train.u-data.train.u(1)  ,data.train.tVec)+data.train.u(1);

    dataFilt.valid.y   = lsim(LPfilt,data.valid.y-data.valid.y(1)    ,data.valid.tVec)  +data.valid.y(1);
    dataFilt.valid.u    = lsim(LPfilt,data.valid.u-data.valid.u(1)      ,data.valid.tVec)  +data.valid.u(1);

    dataFilt.full.y  = lsim(LPfilt,data.full.y-data.full.y(1)  ,data.full.tVec) +data.full.y(1);
    dataFilt.full.u   = lsim(LPfilt,data.full.u-data.full.u(1)    ,data.full.tVec) +data.full.u(1);
    
%% Make iddata sets from the data
    data.train.id = iddata(data.train.y,data.train.u,Ts);
    data.valid.id = iddata(data.valid.y,data.valid.u,Ts);
    data.trans.id = iddata(data.trans.y,data.trans.u,Ts);
    data.full.id  = iddata(data.full.y, data.full.u, Ts);
    
    dataFilt.train.id = iddata(dataFilt.train.y,dataFilt.train.u,Ts);
    dataFilt.valid.id = iddata(dataFilt.valid.y,dataFilt.valid.u,Ts);
    dataFilt.trans.id = iddata(dataFilt.trans.y,dataFilt.trans.u,Ts);
    dataFilt.full.id  = iddata(dataFilt.full.y, dataFilt.full.u, Ts);



end




























