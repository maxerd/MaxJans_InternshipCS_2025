function u = chooseInput(tVec,input,Watt,varargin)
% u = chooseInput(tVec,input,Watt,varargin)
% 
% Input vector generation with 5 rows
% 
%% Inputs
%   tVec     --> Time vector related to the input vector
%   input    --> Input name, always in combination with Watt, options:
%                 -> HT1 - Heater 1
%                 -> HT2 - Heater 2
%                 -> HT3 - Heater 3
%                 -> HT4 - Heater 4
%                 -> HT5 - Heater 5
%   Watt     --> Input value given to the input, can be constant or time varying
%                       Always in combination with input
%   varargin --> Space for additional input definitions,
%                       Always in pairs, inputName and wattage
%                 -> input
%                 -> Watt
%
%% Outputs
%   u --> Generated input vector
%
%% Example
%    Heater 3 is used with a constant 4W
%      u = chooseInput(tVec,'HT3',4); 
%    
%    Heater 1, 3 and 5 are used with a constant 4W
%      u = chooseInput(tVec,'HT1',4,'HT3',4,'HT5',4); 
%    
%    All heater are used with a constant 4W
%      u = chooseInput(tVec,'all',4); 
%    
%    None of the heaters are active
%      u = chooseInput(tVec,'none',0); 
%    
%    Heater 3 is active with a sine wave input
%      u = chooseInput(tVec,'HT3',2*sin(tVec)+2); 
%
%% Made by Max Jans, TU/e, File creation 13-08-2025
%                          Last updated  13-08-2025
% 

    
    if and(size(tVec,2)~=size(Watt,2),length(Watt)~=1)
        error('Input vector sizes to not match, please correct')
    end

    inputName{1}   = input;
    inputWatt(1,:) = Watt;
    nInputs = 1;
    if ~isempty(varargin)
        if ~rem(size(varargin,2),2)
            nInputs = size(varargin,2)/2+1;
            for i=1:nInputs-1
                if and(size(tVec,2)~=size(varargin{2*(i)},2),length(Watt)~=1)
                    error('Input vector sizes to not match, please correct')
                end
                inputName{i+1} = varargin{2*(i)-1};
                inputWatt(i+1,:) = varargin{2*(i)};
            end
        else
            warning('Not enough inputs given, only first input will be considered')
        end
    end
    u = zeros(size(tVec,2),5);
    for i=1:nInputs
        if strcmp(inputName{i},'HT1')
            u(:,1) = inputWatt(i,:).*ones(size(tVec,2),1)';
        elseif strcmp(inputName{i},'HT2')
            u(:,2) = inputWatt(i,:).*ones(size(tVec,2),1)';
        elseif strcmp(inputName{i},'HT3')
            u(:,3) = inputWatt(i,:).*ones(size(tVec,2),1)';
        elseif strcmp(inputName{i},'HT4')
            u(:,4) = inputWatt(i,:).*ones(size(tVec,2),1)';
        elseif strcmp(inputName{i},'HT5')
            u(:,5) = inputWatt(i,:).*ones(size(tVec,2),1)';
        elseif strcmp(inputName{i},'none')
            warning('The "none" option discards the "Watt" input, act accordingly')
            u      = zeros(size(tVec,2),5)';
        elseif strcmp(inputName{i},'all')
            warning('The "all" option only works reliably for time-invariant watt inputs, act accordingly')
            u      = inputWatt(i,:).*ones(size(tVec,2),5)';
        else
            warning('At least one of the input names is unknown')
        end
    end

    u = u';

end