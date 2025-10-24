function [data, dataFilt] = makeOLdata(G,tMeas,dist,v,Ts)

%% Some user inputs
makeValSet  = true;  % [true/false] Is a seperate validation set needed?
removeTrans = false; % [true/false] Does the transient need to be removed

%% Define some simulation variables
fs = 1/Ts;     % [Hz]    Sampling frequency

T_amb = 23; % [degC] Ambient Temperature
T0 = 23; % [degC] Inital Temperature

N = tMeas*3600*fs; % [-] Measurement samples
tVec = linspace(fs,tMeas*3600,N); % [s] Time vector

%% Define the noise on the signal
% v     = 0.15.*randn(N,1);
v_val = v;

%% Simulate the systems

disp('Running simulation')

% Ambient temperature vector, constant ambient temperature for simulations
    ambVec = T_amb.*ones(1,length(dist));

% Define the initial temperatures (states) for the simulation, also constant.
    x0_OL = T0.*ones(size(G.A,1),1);

% Run the open loop simulation, with the previously defined noise signal added
    y_OL     = lsim(G,[dist    ;ambVec],tVec,x0_OL)+v;
    y_OL_val = lsim(G,[dist;ambVec],tVec,x0_OL)+v_val;

% First order parametric approximation, using time data
    sysID.par.firstAprrox.OL.raw = step_sysID(dist,zeros(size(dist)),y_OL-ambVec',tVec,10000*fs);

% Remove the transient using the first order approximation
    y_approx1     = lsim(sysID.par.firstAprrox.OL.raw,[dist],tVec);
    y_approx1_val = lsim(sysID.par.firstAprrox.OL.raw,[dist],tVec);
    % figure;plot(tVec,y_OL-y_approx1)
    % figure;plot(tVec,y_OL,tVec,y_approx1)
    % y_OL     = ((y_OL-y_approx1));
    % figure;plot(tVec,y_OL);grid minor
    % y_OL_val = (y_OL_val-y_approx1_val);

%% Define the in and output vector sizes

% Define the range of the data that is used for identification and validation
if removeTrans
    idx1_OL = 8*3600*fs;
else
    idx1_OL = 1;
end
if makeValSet
    idx2_OL = N-(N-6*3600*fs)*0.2;
else
    idx2_OL = N;
end

idxRange_OL_trns = 1:idx1_OL+1;
idxRange_OL      = idx1_OL:idx2_OL;
idxRange_OL_val  = (idx2_OL-1):N;

%% Define the in and output vectors, for easy and consistent use
% Open loop
data.trans.y = y_OL(idxRange_OL_trns)-ambVec(idxRange_OL_trns)';
data.trans.u = dist(idxRange_OL_trns)';

data.train.y = y_OL(idxRange_OL)-ambVec(idxRange_OL)';
data.train.u = dist(idxRange_OL)';

data.valid.y = y_OL_val(idxRange_OL_val)-ambVec(idxRange_OL_val)';
data.valid.u = dist(idxRange_OL_val)';

data.full.y  = y_OL-ambVec';
data.full.u  = dist';

% Define the time vectors
data.trans.tVec = tVec(idxRange_OL_trns);
data.train.tVec = tVec(idxRange_OL);
data.valid.tVec = tVec(idxRange_OL_val);
data.full.tVec  = tVec;

%% Filter the data vectors to remove some of the effect of noise
% High order lowpass filter, to filter away all the noise above the 
% multisine exitation frequency
    LPfilt = orderLP(0.01,2,0,0,0.7)^2;

% Actually filter the data
    % Open loop data
    dataFilt.trans.y = lsim(LPfilt,data.trans.y-data.trans.y(1),data.trans.tVec)  +data.trans.y(1);
    dataFilt.trans.u = lsim(LPfilt,data.trans.u-data.trans.u(1),data.trans.tVec)  +data.trans.u(1);

    dataFilt.train.y = lsim(LPfilt,data.train.y-data.train.y(1),data.train.tVec)  +data.train.y(1);
    dataFilt.train.u = lsim(LPfilt,data.train.u-data.train.u(1),data.train.tVec)  +data.train.u(1);

    dataFilt.valid.y   = lsim(LPfilt,data.valid.y-data.valid.y(1),data.valid.tVec)+data.valid.y(1);
    dataFilt.valid.u   = lsim(LPfilt,data.valid.u-data.valid.u(1),data.valid.tVec)+data.valid.u(1);

    dataFilt.full.y  = lsim(LPfilt,data.full.y-data.full.y(1),data.full.tVec)     +data.full.y(1);
    dataFilt.full.u  = lsim(LPfilt,data.full.u-data.full.u(1),data.full.tVec)     +data.full.u(1);

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


