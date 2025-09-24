function [H,coherence] = makeClosedLoopFRF(output,disturbance,C,fs)

    % Ts = 1/4000;
    % fs = 1/Ts;
    nfft = ceil(length(output)/10);
    window = hann(nfft);
    noverlap = nfft/2;
    
    res = fs/nfft;
    disp("Estimation resolution: "+ num2str(res)+"Hz")
    
    [psd_Noise,hz_Noise] = cpsd(disturbance,disturbance,window,noverlap,nfft,fs);
    [cpsd_HoutNoise,hz_HoutNoise] = cpsd(output,disturbance,window,noverlap,nfft,fs);
    Sud = frd(cpsd_HoutNoise,hz_HoutNoise*2*pi);
    Sdd = frd(psd_Noise,hz_Noise*2*pi);
    
    S = Sud/Sdd;
    H = (C^-1)*((S^-1)-1);

    [coh,hz] = mscohere(disturbance,output,window,noverlap,nfft,fs);

    coherence = frd(coh,hz*2*pi);
end