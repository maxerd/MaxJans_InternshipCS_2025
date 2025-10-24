
main_lumpedSystem_water_V3

%%
clear mag phase vWat magTemp phaseTemp wOutTemp
% vWat = linspace(0,0.005,20);
vWat = linspace(0,1,200);


for i=1:length(vWat)
    disp(['Iteration: ',num2str(i)])
    sysTMC_sample{i} = usubs(sysTMC(216,3),'waterVel',vWat(i));
    % sysTMC_sample{i} = usample(sysTMC(218,3));
    [magTemp,phaseTemp,wOutTemp] = bode(sysTMC_sample{i},logspace(-6,-1,100));
    mag(:,i) = squeeze(magTemp);
    phase(:,i) = squeeze(phaseTemp);
end

%%


figure(101);
    subplot(211)
        surf(vWat,wOutTemp,db(mag))
        % shading flat
        set(gca, 'YScale', 'log')
    subplot(212)
        surf(vWat,wOutTemp,phase)
        set(gca, 'YScale', 'log')







