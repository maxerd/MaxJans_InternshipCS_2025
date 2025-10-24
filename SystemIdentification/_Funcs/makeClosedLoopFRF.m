function [H,coherence] = makeClosedLoopFRF(output,disturbance,C,fs)
% [H,coherence] = makeClosedLoopFRF(output,disturbance,C,fs)
% 
% Make the Frequency Reponse Function from open loop (direct) data
% 
%% Inputs
%   output      --> The output for the identification
%                       -> For this method either,
%                           - Error signal
%                           - Plant input (U)
%   disturbance --> The input for the identification (Exitation signal)
%                       -> For this method either,
%                           - Reference signal
%                           - Plant input disturbance (d)
%   C           --> Controller that was used during the identfication
%   fs          --> The sampling frequency of the data
%
%% Outputs
%   H           --> The resulting FRF
%   coherence   --> The coherence function correponding to the generated FRF
%
%% Information
%  This function uses the sensitivity function to estimate the FRF,
%  therefore both the transfer from y->e and d->u can be used. 
%       (Classical Indirect identification)
% 
%  The coherence function can be deceiving for closed loop FRFs, dont put
%  too much emphasis on it
%
%% Made by Max Jans, TU/e, File creation xx-09-2024
%                         Last updated  15-10-2025

    % Define how many windows are used for the identification process
    nWindow = 5;

    % Calculate how many samples are in each window
    nfft = ceil(length(output)/nWindow);

    % Define the window and the window overlap
    window = hann(nfft);
    noverlap = nfft/2;
    
    % Give the user the FRF resolution
    res = fs/nfft;
    disp("Estimation resolution: "+ num2str(res)+"Hz")
    
    % Get the PS (power spectrum) of the input signal
    [psd_Noise,hz_Noise] = cpsd(disturbance,disturbance,window,noverlap,nfft,fs);

    % Get the CPS (cross power spectrum) of the input to the output
    [cpsd_HoutNoise,hz_HoutNoise] = cpsd(output,disturbance,window,noverlap,nfft,fs);

    % Define the power spectra as frequency response data
    Sud = frd(cpsd_HoutNoise,hz_HoutNoise*2*pi);
    Sdd = frd(psd_Noise,hz_Noise*2*pi);
    
    % Get the FRF model from the power spectra
    S = Sud/Sdd;
    H = (C^-1)*((S^-1)-1);

    % Make the coherence, also in a frequency reponse data format
    [coh,hz] = mscohere(disturbance,output,window,noverlap,nfft,fs);
    coherence = frd(coh,hz*2*pi);
end