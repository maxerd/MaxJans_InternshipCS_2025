function [G_poly,T_poly] = sysID_LPM(outputData,inputData,fs,n,R,varargin)
%% [G_poly,T_poly] = sysID_LPM(outputData,inputData,fs,n,R,varargin)
%
% Perform system identification using a local polynomial method (LPM)
%       For SISO system only
% 
%% Inputs
%   outputData  --> Plant output data
%   inputData   --> Plant input data
%   fs          --> Sampling frequency of input/output data
%   n           --> Amount of samples forward and backwards that the fit is made on
%   R           --> Order of the polynomial for the fit
%   varargin    --> Option for period data
%                   -> Use a number for the amount of periods
%
%% Outputs
%   G_poly --> Identified plant using the LPM
%   T_poly --> Identified transient using the LPM
%
%% Example
%   For identification of a system with random input, N=20 input is
%   repeated 5 times in measurement
%       [G_poly,T_poly] = sysID_LPM(outputData,randn(100,1),fs,n,R,5)
%
%   !!!Not yet implemented!!!
%     For identification of a system with random input, N=20 input is
%     repeated 5 times in measurement
%         [G_poly,T_poly] = sysID_LPM(outputData_20x5,randn(50,5),fs,n,R)
%   !!!Not yet implemented!!!
%
%% References
%  [3] Evers, E. (2020). Identification and Active Thermomechanical Control in Precision Mechatronics.
%  [2] Pintelon-system-identification-a-frequency-domain-approach-2ed
%  [1] J. Schoukens, R. Pintelon, Y. Rolain - Identification of linear dynamic systems
%
%% Made by Max Jans, TU/e, File creation 18-08-2025
%                          Last updated  18-08-2025

%% Quick check
if ~((2*n+1)>=(R+1)*(1+1))
    error('Results will most likely be useless, therefore script is stopped')
end


%% Check if data is in correct format
if size(inputData,1)~=1
    disp('Input data is going to be transposed')
    inputData  = inputData';
    if size(outputData,1)~=1
    disp('Output data is going to be transposed')
        outputData = outputData';
    end

elseif size(inputData,2)~=1
    disp('Data is fine!')
    if size(outputData,2)~=1
        outputData = outputData';
    else
        disp('Output data is going to be transposed')
    end
else
    % For future: simply detect the data shape and assume periodicity
    warning('Data needs to be reshaped, please check answers carefully.')
    inputData  = reshape(inputData,1,[]);
    outputData  = reshape(outputData,1,[]);
end

u = inputData;
y = outputData;

%% Check for periodicity in the data
periods = 1;
if ~isempty(varargin)
    disp(['Periods in in/output data: ', num2str(varargin{1})])
    periods = varargin{1};
    u = mean(reshape(inputData ,[],periods),2);
    y = mean(reshape(outputData,[],periods),2);
end

%% System identification

% Define some needed variables
N = length(u);
r = [-n:n];
k = [n:N/2-n]+1;

% DFT of the in/output data
U = fft(u);
Y = fft(y');

% Perform least squares estimation of the needed parameters [2]
% Kn = zeros((R+1),(2*n+1));
for i=k
    for j=1:(2*n+1)
        % for l = 1:(R+1)
        %     Kr(l,j) = r(j)^(l-1);%
        % end
        Kr(:,j) = [1; r(j); r(j)^2; r(j)^3; r(j)^4];
        Kn(:,j) = [kron(Kr(:,j),U(i+r(j))); Kr(:,j)];
    end
    % size(Y(i+r))
    % size(Kn)
    par(i-n,:) = Y(i+r)*Kn'*inv((Kn*Kn'));
end

%% Ensemble system
% fVec = (k*fs/N);
fVec = linspace(fs/N,fs/2,N/2);
fVec = fVec(k)*periods;

% Make FRD model from data for easy use with 'bode()' etc
G_poly = frd(par(:,1),fVec*2*pi);
T_poly = frd(par(:,5),fVec*2*pi);



end