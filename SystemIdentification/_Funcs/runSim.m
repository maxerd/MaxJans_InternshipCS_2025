try 
    sysTMC;
    disp('Model found')
catch
    warning('Model not found, running model generation')
    main_lumpedSystem_water_V3
    warning('Model was not found, model generation done')
end

%%

Ts_sim = data.tVec(2)-data.tVec(1);
Tend = data.tVec(end);

simDat.y_e  = timetable(seconds(data.tVec)',data.y_e' );
simDat.amb  = timetable(seconds(data.tVec)',data.amb' );
simDat.d    = timetable(seconds(data.tVec)',data.d'   );
simDat.useC = timetable(seconds(data.tVec)',data.useC');
simDat.ref  = timetable(seconds(data.tVec)',data.ref' );

out = sim('./allModel_V1slx.slx');


% Some user inputs
makeValSet  = true;  % [true/false] Is a seperate validation set needed?
removeTrans = true; % [true/false] Does the transient need to be removed

%% Define the in and output vector sizes

% Define the range of the data that is used for identification and validation
transTime = 3.5; % [h] transient time
valSplit = 0.2; % How much of the data is validation set

if removeTrans
    idx1_CL = transTime*3600*fs;
else
    idx1_CL = 1;
end
if makeValSet
    idx2_CL = N-(N-transTime*3600*fs)*valSplit;
else
    idx2_CL = N;
end

idxRange_CL_trns = 1:idx1_CL+1;
idxRange_CL      = idx1_CL:idx2_CL;
idxRange_CL_val  = (idx2_CL):N;

%%
dataOut.y      = out.y.Data;
dataOut.amb    = out.amb.Data;
dataOut.y_real = out.y_real.Data;
dataOut.u      = out.u.Data;
dataOut.dist   = out.dist.Data;
dataOut.Cout   = out.Cout.Data;
dataOut.err    = out.err.Data;
dataOut.ref    = out.ref.Data;
dataOut.tVec   = out.y.Time;

%% Define the in and output vectors, for easy and consistent use
sysID.data.CL.trans.y = (out.y.Data(idxRange_CL_trns)-out.amb.Data(idxRange_CL_trns));
sysID.data.CL.trans.u = (out.u.Data(idxRange_CL_trns));
sysID.data.CL.trans.d = (out.dist.Data(idxRange_CL_trns));
sysID.data.CL.trans.e = (out.err.Data(idxRange_CL_trns));
sysID.data.CL.trans.r = (out.ref.Data(idxRange_CL_trns));

sysID.data.CL.train.y = (out.y.Data(idxRange_CL)-out.amb.Data(idxRange_CL));
sysID.data.CL.train.u = (out.u.Data(idxRange_CL));
sysID.data.CL.train.d = (out.dist.Data(idxRange_CL));
sysID.data.CL.train.e = (out.err.Data(idxRange_CL));
sysID.data.CL.train.r = (out.ref.Data(idxRange_CL));

sysID.data.CL.valid.y = (out.y.Data(idxRange_CL_val)-out.amb.Data(idxRange_CL_val));
sysID.data.CL.valid.u = (out.u.Data(idxRange_CL_val));
sysID.data.CL.valid.d = (out.dist.Data(idxRange_CL_val));
sysID.data.CL.valid.e = (out.err.Data(idxRange_CL_val));
sysID.data.CL.valid.r = (out.ref.Data(idxRange_CL_val));

sysID.data.CL.full.y = (out.y.Data-out.amb.Data);
sysID.data.CL.full.u = (out.u.Data);
sysID.data.CL.full.d = (out.dist.Data);
sysID.data.CL.full.e = (out.err.Data);
sysID.data.CL.full.r = (out.ref.Data);

% Define the time vectors
sysID.data.CL.trans.tVec = out.y.Time(idxRange_CL_trns);
sysID.data.CL.train.tVec = out.y.Time(idxRange_CL);
sysID.data.CL.valid.tVec = out.y.Time(idxRange_CL_val);
sysID.data.CL.full.tVec  = out.y.Time;















































