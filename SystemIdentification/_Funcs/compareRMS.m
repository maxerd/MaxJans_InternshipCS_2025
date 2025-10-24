function varargout = compareRMS(data,dataNames,sys,varargin)
% compareRMS(data,dataNames,sys,varargin)
% 
% Make a plot such as with Matlab's own 'compare.m', but with a
% (non-normalized) RMS metric
% 
%% Inputs
%   data      --> The iddata that the comparison needs to be performed on
%   dataNames --> The data names that should be displayed in the legend
%                 -> The 'data' legend entity is added by default, so the
%                    amount of names should be the same as the amount of systems
%   sys       --> The systems that need to be compared
%                 -> Can be 'ss', 'tf', 'frd' and 'idss'
%   varargin  --> More systems and kstep
%                 -> More systems like in the input 'sys' can be added
%                 -> kstep can be added -> The amount of prediction steps
%                    -> simulation   = inf
%                    -> x-step ahead = x
%
%% Example
% dataSet = iddata(y,u,Ts)
% 
% Compare the simulation of two models, from 'Method 1' and 'Method 2'
%    figure(1);clf;
%       compareRMS(dataSet,{'Method 1','Method 2'},mdl_1,mdl_2,inf);
%           title('Comparison of Method 1 and Method 2')
%    
% Compare the 10-step ahead prediction of two models, from 'Method 1' and 'Method 2'
%    figure(1);clf;
%       compareRMS(dataSet,{'Method 1','Method 2'},mdl_1,mdl_2,10);
%           title('Comparison of Method 1 and Method 2')
%    
% Compare the 10-step ahead prediction of three models
%    figure(1);clf;
%       compareRMS(dataSet,{'Method 1','Method 2','Method 3'},mdl_1,mdl_2,mdl_3,10);
%           title('Comparison of Method 1, Method 2 and Method 3')
%    
%% Made by Max Jans, TU/e, File creation 26-09-2025
%                         Last updated  26-09-2025
% 
% See also compare, iddata, rmse

%% Input handling
    % Default steps ahead --> Simulation
    xStep = inf;
    
    % There should always be one system present, as defined in the inputs
    systems{1} = sys;
    for j=1:(nargin-3)
    
        % Define the systems that need to be compared
        if or(or(strcmp(class(varargin{j}),'ss'),strcmp(class(varargin{j}),'tf')),or(strcmp(class(varargin{j}),'frd'),strcmp(class(varargin{j}),'idss')))
            systems{j+1} = varargin{j};
        end
    
        % Define the amount of steps that need to be predicted
        if strcmp(class(varargin{j}),'double')
            xStep = varargin{j};
        end
    end
    
%% Use the 'compare.m' function to get the prediction of the data with the right inital conditions
    for k=1:length(systems)
        [ymod{k},~,~] = compare(data,systems{k},xStep);
        compData(:,k) = ymod{k}.y;
    end
    
    % Define the time vector from the iddata structure
    timeVec = 0:data.Ts:(length(data.y)*data.Ts-data.Ts);
    
    % The standard matlab colors to give to the plots
    mlColors = orderedcolors("gem");
    
%% Optional output definition, also removing the plotting if true
    plotOpt = 1;
    if nargout > 0
        plotOpt = 0;
        varargout{1} = compData;
    end
    
    
%% If needed, plot the data
    if plotOpt
        % First plot of the given data, in grey
        p = plot(timeVec,data.y);hold on;grid minor
        p(1).Color = [0.5 0.5 0.5];
        
        % Initialize the legend data
        leg = {'Data'};
        
        % Plot the prediction for every system
        for i=1:size(ymod,2)
            p = plot(timeVec,(ymod{i}.y));
        
            % Define the color that this loop's plot should have, restarting at the
            % first one when the last one is reached, just like native matlab functions
                colorRGB = mlColors(rem(i-1,7)+1,:);
                p.Color = colorRGB;
        
            % Option to add the gap metric to the legend instead of the RMSE, this
            % does decrease performance, and the systems are tested against the
            % first input system
                % [gap,~] = gapmetric(systems{1},systems{i});
                % leg{i+1} = [dataNames{i},': ',num2str(gap,3),' RMSE'];
        
            % Make the legend from the dataNames structure
                leg{i+1} = [dataNames{i},': ',num2str(rmse(ymod{i}.y,data.y),3),' RMSE'];
        end
        
        % Standard labels for current project, can be changed after use of the
        % function
            xlabel('Time [s]')
            ylabel('Temperature [degC]')
            legend(leg)
    end

end




