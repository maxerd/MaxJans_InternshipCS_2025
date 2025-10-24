function [H,coherence] = makeOpenLoopFRF(output,disturbance,fs)

    Ts = 1/fs;
    nfft = ceil(length(output)/5);
    window = hann(nfft);
    noverlap = nfft/2;
    
    res = fs/nfft;
    disp("Estimation resolution: "+ num2str(res)+"Hz")
    
    [psd_dis,hz_dis] = cpsd(disturbance,disturbance,window,noverlap,nfft,fs);
    [cpsd_outDis,hz_outDis] = cpsd(output,disturbance,window,noverlap,nfft,fs);
    Syd = frd(cpsd_outDis,hz_outDis*2*pi);
    Sdd = frd(psd_dis,hz_dis*2*pi);
    
    H = Syd/Sdd;

    [coh,hz] = mscohere(disturbance,output,window,noverlap,nfft,fs);

    coherence = frd(coh,hz);
end