clc

data.tVec = 0:0.1:20000;
data.y_e  = 0.15.*randn(size(data.tVec));
data.d    = 1.*randn(size(data.tVec))+20;
data.ref  = 40.8.*ones(size(data.tVec));

data = chooseLoop(data,'partial',3);

runSim

figure(7);clf;plot(dataOut.tVec,dataOut.y)

















