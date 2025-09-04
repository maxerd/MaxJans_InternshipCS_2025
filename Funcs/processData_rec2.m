function data = processData_rec2(path,opts)
% data = processData_rec1(path,opts)
% 
% Process measurement data for the TU/e, CS-group, Thermal Mass Control Setup
% 
%% Inputs
%   path --> The absolute path of the raw data file (.mat file)
%   opts --> Processing options (struct)
%            -> opts.save,     boolean, to save or not to save?
%                                          if false, output is a struct
%                                          if true, no direct output is a given
%            -> opts.savePath, string,  filepath+filename to save the data to
%            -> opts.downSamp, double,  downsample by taking only the n-th data points
%
%% Outputs
%   data --> Processed data
%
%% Example
%    Still needs to be added
%      opts.save       = false;
%      data = processData_rec2('D:\_data\NoiseMeasurement.mat',opts);
%
%% Made by Max Jans, TU/e, File creation 19-08-2025
%                          Last updated  19-08-2025

%% Load in data
try
    rawData = load(path);
catch
    error('Error in loading data, check file name and extension (need .mat)')
end

%% Split data
% Should also add a trigger possibility to each signal?

try 
    NdownSample = opts.downSample;
catch
end

dataNameAll = split(path,'\');
dataName = dataNameAll{end}(1:end-4);

%% Get the index list of the signals
% This is used to automatically select the index of the named signals
%   The signals change index if the "Recorder" in dSpace is (slightly)
%   changed, making data processing difficult otherwise

nameList = {'T',...
            'Temp[0]',...
            'Temp[1]',...
            'Temp[2]',...
            'Flow',...
            'ID',...
            'Recv Frames',...
            'Status',...
            'Temp flow',...
            'Temp tank',...
            'Timestamp high',...
            'Timestamp low',...
            'In1',...
            'Out1',...
            'Out1[0]',...
            'T_SP'};

idxList = zeros(length(nameList),1);
for j=1:length(eval(['rawData.' dataName '.Y']))
    for i=1:length(nameList)
        if strcmp(eval(['rawData.' dataName '.Y(' num2str(j)    ').Name']),nameList(i))
            idxList(i) = j;
        end
    end
end


%% Get the signal index connected to the right signal names again

idx.periods = idxList(13);

idx.Watt    = idxList(14);
idx.tempHT  = idxList(3);
idx.tempAmb = idxList(4);
idx.tempTM  = idxList(2);
idx.tempTMf = idxList(1);
idx.flow    = idxList(5);
idx.tempFlw = idxList(9);
idx.ID      = idxList(6);
idx.rcvFrms = idxList(7);
idx.status  = idxList(8);
idx.tempTnk = idxList(10);
idx.error   = idxList(15);
idx.SP      = idxList(16);

% The "Watt" input has the tendency to dissapear from the recorder, this if
% statement is there to avoid errors
if idx.Watt
    namedData.Watt = eval(['rawData.' dataName '.Y(' num2str(idx.Watt)    ').Data']);
else
    namedData.Watt = nan;
    warning('No "Watt" signal detected, please act accordingly')
end

namedData.tVec    = eval(['rawData.' dataName '.X.Data']);
namedData.tempHT  = eval(['rawData.' dataName '.Y(' num2str(idx.tempHT)  ').Data']);
namedData.tempAmb = eval(['rawData.' dataName '.Y(' num2str(idx.tempAmb) ').Data']);
namedData.tempTM  = eval(['rawData.' dataName '.Y(' num2str(idx.tempTM)  ').Data']);
if idx.tempTMf
    namedData.tempTMf = eval(['rawData.' dataName '.Y(' num2str(idx.tempTMf) ').Data']);
else
    namedData.tempTMf = nan;
    warning('No "filtered" TM data')
end
namedData.flow    = eval(['rawData.' dataName '.Y(' num2str(idx.flow)    ').Data']);
namedData.tempFlw = eval(['rawData.' dataName '.Y(' num2str(idx.tempFlw) ').Data']);
namedData.ID      = eval(['rawData.' dataName '.Y(' num2str(idx.ID)      ').Data']);
namedData.rcvFrms = eval(['rawData.' dataName '.Y(' num2str(idx.rcvFrms) ').Data']);
namedData.status  = eval(['rawData.' dataName '.Y(' num2str(idx.status)  ').Data']);
namedData.tempTnk = eval(['rawData.' dataName '.Y(' num2str(idx.tempTnk) ').Data']);

if idx.tempTMf
    namedData.error = eval(['rawData.' dataName '.Y(' num2str(idx.error) ').Data']);
    namedData.SP    = eval(['rawData.' dataName '.Y(' num2str(idx.SP)    ').Data']);
else
    namedData.error   = nan;
    namedData.SP      = nan;
    warning('No setpoint data')
end

if opts.save
    if ~exist('NdownSample')
        tVec    = namedData.tVec;
        Watt    = namedData.Watt;
        tempHT  = namedData.tempHT;
        tempAmb = namedData.tempAmb;
        tempTM  = namedData.tempTM;
        tempTMf = namedData.tempTMf;
        error   = namedData.error;
        SP      = namedData.SP;

        Ts      = namedData.tVec(3)-namedData.tVec(2);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Data coming from raspberry pi, if error occurs here check UDP connection
        flow    = namedData.flow   ;
        tempFlw = namedData.tempFlw;
        ID      = namedData.ID     ;
        rcvFrms = namedData.rcvFrms;
        status  = namedData.status ;
        tempTnk = namedData.tempTnk;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    else
        tVec    = (mean(reshape(namedData.tVec(1:end-rem(length(namedData.tVec),opts.downSample)),      opts.downSample,[])));
        Watt    = (mean(reshape(namedData.Watt(1:end-rem(length(namedData.Watt),opts.downSample)),      opts.downSample,[])));
        tempHT  = (mean(reshape(namedData.tempHT(1:end-rem(length(namedData.tempHT),opts.downSample)),  opts.downSample,[])));
        tempAmb = (mean(reshape(namedData.tempAmb(1:end-rem(length(namedData.tempAmb),opts.downSample)),opts.downSample,[])));
        tempTM  = (mean(reshape(namedData.tempTM(1:end-rem(length(namedData.tempTM),opts.downSample)),  opts.downSample,[])));
        tempTMf = (mean(reshape(namedData.tempTMf(1:end-rem(length(namedData.tempTMf),opts.downSample)),opts.downSample,[])));
        error   = (mean(reshape(namedData.error(1:end-rem(length(namedData.error),opts.downSample)),opts.downSample,[])));
        SP   = (mean(reshape(namedData.SP(1:end-rem(length(namedData.SP),opts.downSample)),opts.downSample,[])));

        Ts      = (namedData.tVec(3)-namedData.tVec(2))*opts.downSample;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Data coming from raspberry pi, if error occurs here check UDP connection
        flow    = (mean(reshape(namedData.flow(1:end-rem(length(namedData.flow),opts.downSample)),      opts.downSample,[])));
        tempFlw = (mean(reshape(namedData.tempFlw(1:end-rem(length(namedData.tempFlw),opts.downSample)),opts.downSample,[])));
        ID      = (mean(reshape(namedData.ID(1:end-rem(length(namedData.ID),opts.downSample)),          opts.downSample,[])));
        rcvFrms = (mean(reshape(namedData.rcvFrms(1:end-rem(length(namedData.rcvFrms),opts.downSample)),opts.downSample,[])));
        status  = (mean(reshape(namedData.status(1:end-rem(length(namedData.status),opts.downSample)),  opts.downSample,[])));
        tempTnk = (mean(reshape(namedData.tempTnk(1:end-rem(length(namedData.tempTnk),opts.downSample)),opts.downSample,[])));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
else
    if ~exist('NdownSample')
        data.tVec    = namedData.tVec;
        data.Watt    = namedData.Watt;
        data.tempHT  = namedData.tempHT;
        data.tempAmb = namedData.tempAmb;
        data.tempTM  = namedData.tempTM;
        data.tempTMf = namedData.tempTMf;
        data.error   = namedData.error;
        data.SP      = namedData.SP;

        data.Ts      = namedData.tVec(3)-namedData.tVec(2);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Data coming from raspberry pi, if error occurs here check UDP connection
        data.flow    = namedData.flow   ;
        data.tempFlw = namedData.tempFlw;
        data.ID      = namedData.ID     ;
        data.rcvFrms = namedData.rcvFrms;
        data.status  = namedData.status ;
        data.tempTnk = namedData.tempTnk;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
        data.tVec    = (mean(reshape(namedData.tVec(1:end-rem(length(namedData.tVec),opts.downSample)),      opts.downSample,[])));
        data.Watt    = (mean(reshape(namedData.Watt(1:end-rem(length(namedData.Watt),opts.downSample)),      opts.downSample,[])));
        data.tempHT  = (mean(reshape(namedData.tempHT(1:end-rem(length(namedData.tempHT),opts.downSample)),  opts.downSample,[])));
        data.tempAmb = (mean(reshape(namedData.tempAmb(1:end-rem(length(namedData.tempAmb),opts.downSample)),opts.downSample,[])));
        data.tempTM  = (mean(reshape(namedData.tempTM(1:end-rem(length(namedData.tempTM),opts.downSample)),  opts.downSample,[])));
        data.tempTMf = (mean(reshape(namedData.tempTMf(1:end-rem(length(namedData.tempTMf),opts.downSample)),opts.downSample,[])));
        data.error   = (mean(reshape(namedData.error(1:end-rem(length(namedData.error),opts.downSample)),opts.downSample,[])));
        data.SP      = (mean(reshape(namedData.SP(1:end-rem(length(namedData.SP),opts.downSample)),opts.downSample,[])));

        data.Ts      = (namedData.tVec(3)-namedData.tVec(2))*opts.downSample;

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % % Data coming from raspberry pi, if error occurs here check UDP connection
        % data.flow    = downsample(namedData.flow   ,opts.downSample);
        % data.tempFlw = downsample(namedData.tempFlw,opts.downSample);
        % data.ID      = downsample(namedData.ID     ,opts.downSample);
        % data.rcvFrms = downsample(namedData.rcvFrms,opts.downSample);
        % data.status  = downsample(namedData.status ,opts.downSample);
        % data.tempTnk = downsample(namedData.tempTnk,opts.downSample);
        data.flow    = (mean(reshape(namedData.flow(1:end-rem(length(namedData.flow),opts.downSample)),      opts.downSample,[])));
        data.tempFlw = (mean(reshape(namedData.tempFlw(1:end-rem(length(namedData.tempFlw),opts.downSample)),opts.downSample,[])));
        data.ID      = (mean(reshape(namedData.ID(1:end-rem(length(namedData.ID),opts.downSample)),          opts.downSample,[])));
        data.rcvFrms = (mean(reshape(namedData.rcvFrms(1:end-rem(length(namedData.rcvFrms),opts.downSample)),opts.downSample,[])));
        data.status  = (mean(reshape(namedData.status(1:end-rem(length(namedData.status),opts.downSample)),  opts.downSample,[])));
        data.tempTnk = (mean(reshape(namedData.tempTnk(1:end-rem(length(namedData.tempTnk),opts.downSample)),opts.downSample,[])));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end

%% Possibility to save the processed data to a file
if opts.save
    % save([opts.savePath,'.mat'],'Ts','tVec','Watt','tempHT','Volt','tempAmb','tempTM','flow','tempFlw','ID','rcvFrms','status','tempTnk')
    save([opts.savePath,'.mat'],'Ts','tVec','Watt','tempHT','tempAmb','tempTM','tempTMf','flow','tempFlw','ID','rcvFrms','status','tempTnk','error','SP')
    data = struct();
    disp('File processed and saved!')
end

end