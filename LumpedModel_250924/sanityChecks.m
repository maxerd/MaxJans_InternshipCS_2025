
clc
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%% Please check these values for correctness %%')
disp('%%   If incorrect, results might be wrong    %%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(' ')

% Total system mass
sumVal = 0;

for i=1:size(TMA,2)
    sumVal = sumVal+TMA(i).mass+TMB(i).mass;
end

disp('%%%%%%%%%%% MASS %%%%%%%%%%%')
disp(['Total mass: ', num2str(sumVal),'g'])
disp(['Real mass: ', num2str(2089.7),'g'])
disp(' ')

% Total convection area
disp('%%%%%%%%%%% AREA %%%%%%%%%%%')
disp(['Total (top) area: ',num2str(sum(LTM.A)*1e6),'mm^2'])
disp(['Approx. (top) area: ',num2str(sum(0.15*0.26)*1.4*1e6),'mm^2'])
disp(' ')

% Conduction
sumVal = 0;

for i=1:size(TMA,2)
    sumVal = sumVal+eval(['LTM' num2str(i) '.GTB']);
end

disp('%%%%%%%%%%% Conduction %%%%%%%%%%%')
disp(['Conduction top-bottom: ',num2str(sumVal),''])
disp(['Approx conduction top-bottom: ',num2str((c.therm.heatCond.alu*0.15*0.26/0.025)),', (+- 10% should be fine, no water-path/holes etc taken into account)'])
disp(' ')
