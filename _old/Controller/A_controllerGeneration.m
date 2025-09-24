
lumpedSystem_V3

%% Select the system that needs to be controlled
% contSys = newSys;
contSys = c2d(sysTMC(32,[3 6]),0.001);

%% Define the system I/O
contSys.InputName =  {'W','Tamb'};
contSys.OutputName = 'T';

%% Generate a PID controller for a user selected bandwidth
Ck = genPID(d2c(contSys(1,1)),0.0025);
% Ck = genPID(d2c(contSys(1,1)),0.005);
Ck.tot.InputName =  'Terr';
Ck.tot.OutputName = 'W';

refSum = sumblk('Terr = Tref - T');

%% Make the closed loop connection
CL = connect(d2c(contSys),Ck.tot,refSum,{'Tref','Tamb'},{'T','W'});

%% Save the controller
% save('PID_01Sep25.mat','Ck','contSys','CL','traj','SF','tVec_C')

%% %% Uncertain water velocity %% %%

lumpedSystem_water_V3

%%
contSys = sysTMC(:,[3 6 7]);
contSys.InputName = {'Watt','Tamb','Twat'};
contSys.OutputName = {'Ttm','Tht','Twatout'};

refSum = sumblk('Terr = Twatref - Twatout');

OLIC = connect((contSys),refSum,{'Twatref','Tamb','Twat','Watt'},{'Terr','Ttm','Tht','Twatout','Watt','Terr','Ttm','Tht'});

[K,CLperf] = musyn(OLIC,3,1);

%%
CLIC = lft(OLIC,K);

tVec = 0:16000;

Twatref = ones( 1,length(tVec) );
Tamb    = zeros(1,length(tVec));
Twat    = zeros(1,length(tVec));

u = [Twatref;Tamb;Twat];

y = lsim(CLIC,u,tVec);

figure;plot(y)







