function u = genMultisine(fs,N,x,maxA,posOnly)
% u = genMultisine(fs,N,x,maxA,posOnly)
% 
% Generate a random-phase mulsine signal which can be used for identification purposes
% 
%% Inputs
%   fs      --> The sampling frequency needed for the signal
%   N       --> The amount of samples in the signal
%   x       --> The percentage of the maximum frequency domain that needs power 
%                   -> If fs/2 is the maximum frequency, the new maximum
%                      will become x*fs/2
%   maxA    --> Maximal amplitude of the signal
%   posOnly --> Give whether the signal needs to be fully positive or not
%                   -> In some applications (such as thermal systems)
%                      negative input can not easily be applied
%
%% Outputs
%   u       --> The generated random-phase multisine signal
%
%% Information
%  Enforcement of the positiveness is achieved by shifting of the
%  signal. 
%       i.e. u_pos = u-min(u);
% 
%  Enforcement of the maximum amplitude is achieved by scaling of the
%  signal
%       i.e. u_maxAmp = u*maxA/max(u);
%
%% Example
%  fs           = 100; % [Hz]  Sampling frequency
%  tMeas        = 8;   % [sec] Measurement time
%  maxAmplitude = 120; % [W]   Maximum exitation signal
%  positiveOnly = 1;   % [-]   Define whether the exitation signal can be negative
%  u            = genMultisine(fs, tMeas*fs, 1, maxAmplitude, positiveOnly);
%
%% Made by Max Jans, TU/e, File creation xx-03-2025
%                         Last updated  15-10-2025

disp('Generating multisine')

% Define the minimum and maximum frequency component in the generated signal
f0 = fs/N;
disp(['Lowest frequency in data: ',num2str(f0),'Hz'])
disp(['Highest frequency in data: ',num2str(x*fs/2),'Hz'])

tVec = 1:N;
A = ones(size(tVec));
phi = randn(N,1)*2*pi;
u = zeros(size(tVec));

maxFreq = floor(x*(N/2-1));

% Actually perform the loop to generate the signal
%   For every time instance, generate a sinusoid with a random phase and
%   add it to the current value.
%       $u(t) = \sum_{j=1}^{\frac{N}{2}-1}A_j\cdot sin(2\pi j\frac{f_0}{f_s}t+\varphi_j)$
%           Where N is the amount of samples in the signal, A_j is the amplitude of
%           each sine wave and \varphi_r the random phase for each wave.
for k=tVec
    for j=1:maxFreq
        u(k) = u(k)+A(j)*sin(2*pi*j*f0/fs*k+phi(j));
    end
end

% Adjusts for positivity and maximum amplitude
if posOnly
    u = u-min(u);
end
u = u*maxA/max(u);

end