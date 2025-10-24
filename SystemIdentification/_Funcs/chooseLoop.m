function data = chooseLoop(data,loop,varargin)
% data = chooseLoop(data,loop,varargin)
% 
% Choose in which mode to run the (linear) simulink model of the TMC-setup
% 
%% Inputs
%   data     --> The data structure defining all the needed simulation variables
%   loop     --> Define in which mode to run the simulation
%                 -> OL - Run the simulation in open loop
%                 -> CL - Run the simulation in closed loop
%                 -> partial - Start the simulation in closed loop, then
%                     switch to open loop, see 'varargin' input
%   varargin --> "p" value, defining what fraction of the simulation is run
%                  in closed loop, only used when loop=='partial'
%                 -> The first 1/p part of the time is done in closed loop
%
%% Example
% data.tVec = 0:0.1:20000;
% data.y_e  = 0.15.*randn(size(data.tVec));
% data.d    = randn(size(data.tVec))+20;
% data.ref  = 40.8.*ones(size(data.tVec));
% 
% data = chooseLoop(data,'partial',2);
%    
% runSim
%    
%% Made by Max Jans, TU/e, File creation 01-10-2025
%                         Last updated  01-10-2025
% 

    N = length(data.tVec);
    if strcmp(loop,'OL')
        % Open loop simulation
        data.useC = zeros(1,floor(N));
        data.ref  = zeros(1,floor(N));
    elseif strcmp(loop,'CL')
        % Closed loop simulation
        data.useC = ones(1,floor(N));
    elseif strcmp(loop,'partial')
        % First 1/p closed loop, after that open loop
        try
            p = varargin{1};
        catch
            warning('No "p" value given, using defaul (p=2)')
            p = 2;
        end
        data.useC = [ones(1,floor(N/p)) zeros(1,ceil(N-N/p))];
    end

end