function [H,coherence] = makeOpenLoopFRF(output,disturbance,fs)
% [H,coherence] = makeOpenLoopFRF(output,disturbance,fs)
% 
% Make the Frequency Reponse Function from open loop (direct) data
% 
%% Inputs
%   output      --> 
%   disturbance --> 
%   fs          --> The sampling frequency of the data
%
%% Outputs
%   H           --> The resulting FRF
%   coherence   --> The coherence function correponding to the generated FRF
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
    [psd_dis,hz_dis] = cpsd(disturbance,disturbance,window,noverlap,nfft,fs);

    % Get the CPS (cross power spectrum) of the input to the output
    [cpsd_outDis,hz_outDis] = cpsd(output,disturbance,window,noverlap,nfft,fs);

    % Define the power spectra as frequency response data
    Syd = frd(cpsd_outDis,hz_outDis*2*pi);
    Sdd = frd(psd_dis,hz_dis*2*pi);
    
    % Get the FRF model from the power spectra
    H = Syd/Sdd;

    % Make the coherence, also in a frequency reponse data format
    [coh,hz] = mscohere(disturbance,output,window,noverlap,nfft,fs);
    coherence = frd(coh,hz);
end




