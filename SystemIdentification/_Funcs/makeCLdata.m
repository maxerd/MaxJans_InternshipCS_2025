function [data, dataFilt, varargout] = makeCLdata(G,tMeas,controllerBW,dist,v,Ts)


%% Some user inputs
makeValSet  = true;  % [true/false] Is a seperate validation set needed?
removeTrans = false; % [true/false] Does the transient need to be removed

%% Adapt the system for use in closed loop simulation
% Define the plant's in and outputnames
    G.inputName = {'HT3','Tamb'};
    G.outputName = 'TM1';

% Generate a controller
    K = genPD(G(1,1),controllerBW);
% Define the plant's in and outputnames
    K.tot.inputName = 'TM1_err';
    K.tot.outputName = 'cOut';

% Sums for interconnection of CL system
    refSum = sumblk('TM1_err = TM1_ref - TM1');
    distSum = sumblk('HT3 = cOut + dist');

% Define the total CL plant model
    P = connect(G,K.tot,refSum,distSum,{'TM1_ref','dist','Tamb'},{'TM1','TM1_err','HT3'});

    if nargin>2
        varargout{1} = K.tot;
    end

%% Define some simulation variables
fs = 1/Ts;     % [Hz]    Sampling frequency

T0 = 23; % [degC] Inital Temperature
Tamb = 23; % [degC] Ambient Temperature

T_ref = 45; % [degC] Reference Temperature

N = round(tMeas*3600*(fs),0); % [-] Measurement samples
% tVec = linspace(fs,tMeas*3600,N); % [s] Time vector
tVec = 0:Ts:(tMeas*3600)-Ts; % [s] Time vector

% Define the noise on the signal
% v     = 0.15.*randn(N,1);
v_val = v;

dist_val = dist;

%% Simulate the systems

disp('Running simulation')

% Ambient temperature vector, constant ambient temperature for simulations
    ambVec = Tamb.*ones(1,length(dist));

% Define the initial temperatures (states) for the simulation, also constant.
    x0_CL = [T0.*ones(size(G.A,1),1); zeros(length(eig(K.tot)),1)];

% The input vector for the CL system consists of the reference temperature,
% the disturbance (/identification) signal and the ambient temperature
    % u     = [T_ref.*ones(1,length(dist))    ;dist    ;ambVec];
    % u_val = [T_ref.*ones(1,length(dist_val));dist_val;ambVec];
    u     = [dist    ;zeros(1,length(dist))    ;ambVec];
    u_val = [dist_val;zeros(1,length(dist_val));ambVec];

% Run the closed loop simulation, again with the previously defined noise signal added
    y_CL     = lsim(P,u    ,tVec,x0_CL)+v;
    y_CL_val = lsim(P,u_val,tVec,x0_CL)+v_val;

%% Define the in and output vector sizes

% Define the range of the data that is used for identification and validation
if removeTrans
    idx1_CL = round(2.5*3600*fs,0);
    % idx1_CL = 20*3600*fs;
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
data.trans.y = y_CL(idxRange_CL_trns,1)-ambVec(idxRange_CL_trns)';
data.trans.u = y_CL(idxRange_CL_trns,3);
data.trans.d = dist(idxRange_CL_trns)';
data.trans.e = y_CL(idxRange_CL_trns,2);
data.trans.r = data.trans.e+data.trans.y;

data.train.y = y_CL(idxRange_CL,1)-ambVec(idxRange_CL)';
data.train.u = y_CL(idxRange_CL,3);
data.train.d = dist(idxRange_CL)';
data.train.e = y_CL(idxRange_CL,2);
data.train.r = data.train.e+data.train.y;

data.valid.y = y_CL(idxRange_CL_val,1)-ambVec(idxRange_CL_val)';
data.valid.u = y_CL(idxRange_CL_val,3);
data.valid.d = dist(idxRange_CL_val)';
data.valid.e = y_CL(idxRange_CL_val,2);
data.valid.r = data.valid.e+data.valid.y;

data.full.y = y_CL(:,1)-ambVec';
data.full.u = y_CL(:,3);
data.full.d = dist';
data.full.e = y_CL(:,2);
data.full.r = data.full.e+data.full.y;

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