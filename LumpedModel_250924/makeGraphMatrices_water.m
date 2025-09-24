
clear connectMap

% Create a Nx by Nz matrix showing the lump numberings
% Example for smaller sizes:
%          [ 3  6  9  12 ]
% numMap = [ 2  5  8  11 ]
%          [ 1  4  7  10 ]
numMap = flip(reshape(ones(1,Nlumps)*(triu(ones(Nlumps))),Nz,Nx));

% Create a Nx by Nz matrix showing the column numberings, used for
%   determing if lumps are touching
% Example for smaller sizes:
%          [ 1  2  3  4 ]
% colMap = [ 1  2  3  4 ]
%          [ 1  2  3  4 ]
colMap = rem(reshape(ones(1,Nlumps)*(triu(ones(Nlumps))),Nx,Nz)-1,9)'+1;

% Initialize the connectMap, which will be translated to the graph's A matrix
%   The entries of the connectMap are the conductances between lumps, however
%   that does mean that to get a series connection of two conductances, we
%   need to do R_tot = 1/G1+1/G2 and G_tot = 1/R_tot
connectMap = zeros(Ny*Nlumps+29); % The 29 lumps are the water/liquid lumps
for i=1:Nlumps
    % Look 7 lumps ahead (neighbors in the x-direction)
    for j=1:(Nlumps-Nz)
        if i==(j+Nz)
            connectMap(i,j) = 1/(1/((TMA(i).GX)*2)+1/((TMA(j).GX)*2));               % For top TM layer in multilayer system
            connectMap(i+Nlumps,j+Nlumps) = 1/(1/((TMB(i).GX)*2)+1/((TMB(j).GX)*2)); % For bottom TM layer in multilayer system
        end
    end
    % Look 1 lump ahead (neighbors in the Z-direction)
    for j=1:(Nlumps-1)
        if and(i==(j+1),colMap(i)==colMap(j))
            connectMap(i,j) = 1/(1/((TMA(i).GZ)*2)+1/((TMA(j).GZ)*2));               % For top TM layer in multilayer system
            connectMap(i+Nlumps,j+Nlumps) = 1/(1/((TMB(i).GZ)*2)+1/((TMB(j).GZ)*2)); % For bottom TM layer in multilayer system
        end
    end
    % Look 1 lump down (neighbors in the Y-direction)
    %   Only from top TM to bottom TM layer
    for j=(Nlumps):2*(Nlumps)
        if i==(j-Nlumps)
            connectMap(i,j) = 1/(1/((TMA(i).GTB)*2)+1/((TMB(j-Nlumps).GTB)*2));      % For multilayer system
        end
    end
    % Look 1 lump down (neighbors in the Y-direction)
    %   Only from bottom TM to heater layer
    for j=2*(Nlumps):3*(Nlumps)
        if i==(j-2*Nlumps)
            connectMap(i+Nlumps,j) = 1/(1/((TMB(i).GTB)*2)+1/((TMC(j-2*Nlumps).GTB)*2)); % For multilayer system
        end
    end
end

waterLumps = [6 13 20:7:48 fliplr([19:7:47]) 18:7:46 fliplr([17:7:45]) 16:7:44 51 58]; % In order of the water flow!
% waterVelocity = 5; %[ml/s]
% waterVelocity = 2; %[ml/s]
% waterVelocity = 1; %[ml/s]
% waterVelocity = 0.5; %[ml/s]
% waterVelocity = 0.1; %[ml/s]
% waterVelocity = 0; %[ml/s]
% waterVelocity = ureal('waterVel',5,'PlusMinus',[-5 5]); %[ml/s]
waterVelocity = waterVelocity*1e-6; % Convert to [m^3/s]
for w = 1:length(waterLumps)
    % 
    if w<=(length(waterLumps)-1)
        % connectMap((end-length(waterLumps)+w),(end-length(waterLumps)+w)+1) = 1/(1/(waterVelocity*c.mech.density.water*c.therm.specHeat.water)+1/(c.therm.heatCond.water*(pi*0.005^4)/0.042));
        connectMap((end-length(waterLumps)+w),(end-length(waterLumps)+w)+1) = c.therm.heatCond.water*(pi*0.005^4)/0.042;
    end
    % if waterVelocity.NominalValue > 0
    if strcmp(class(waterVelocity),'umat')
        if waterVelocity.nominalValue > 0
            ifCondition = 1;
        else
            ifCondition = 0;
        end
    else
        if waterVelocity > 1
            ifCondition = 1;
        else
            ifCondition = 0;
        end
    end
    if ifCondition
        connectMap(waterLumps(w),(end-length(waterLumps)+w))    = 1/(1/(LTM.GWA(waterLumps(w)))+1/(LTM.GWA_cond(waterLumps(w))));
        connectMap(waterLumps(w)+63,(end-length(waterLumps)+w)) = 1/(1/(LTM.GWA(waterLumps(w)+63))+1/(LTM.GWA_cond(waterLumps(w)+63)));
        % connectMap(waterLumps(w),(end-length(waterLumps)+w))    = LTM.GWA(waterLumps(w))+LTM.GWA_cond(waterLumps(w));
        % connectMap(waterLumps(w)+63,(end-length(waterLumps)+w)) = LTM.GWA(waterLumps(w)+63)+LTM.GWA_cond(waterLumps(w)+63);
    else
        connectMap(waterLumps(w),(end-length(waterLumps)+w))    = LTM.GWA_cond(waterLumps(w));
        connectMap(waterLumps(w)+63,(end-length(waterLumps)+w)) = LTM.GWA_cond(waterLumps(w)+63);
    end
end


% Initialize the graph's A matrix
A = zeros(size(connectMap));

% Make the graph's A matrix, including the symmetry
A = connectMap;
A = A'+A;

for w = 1:length(waterLumps)
    % 
    if w<=(length(waterLumps)-1)
        A((end-length(waterLumps)+w)+1,(end-length(waterLumps)+w)) = waterVelocity*c.mech.density.water*c.therm.specHeat.water;
    end
end

% Make the graph's L (laplacian) matrix
L = diag(sum(A,2))-A+diag([LTM.GA zeros(1,(Ny-2)*Nlumps+29)]'+[zeros(1,(Ny)*Nlumps) waterVelocity*c.mech.density.water*c.therm.specHeat.water zeros(1,28)]');


%% Make B matrix for state space model
% No input for the ambient temperature is desired, so last value will
% always to be zero (last state is the ambient temperature state)

% Heater 1 is the left most heater, which is under lump 8-21
%   Since it spans 14 lumps, the input value is divided by 14
% B_heater1 = [zeros(Nz,1);ones(Nz,1)./(2*Nz);ones(Nz,1)./(2*Nz);zeros(Nz,1);zeros(Nz,1);zeros(Nz,1);zeros(Nz,1);zeros(Nz,1);zeros(Nz,1);0];
B_heater1 = [zeros(Nz,1);ones(Nz,1)./(2*Nz);ones(Nz,1)./(2*Nz);zeros(Nz,1);zeros(Nz,1);zeros(Nz,1);zeros(Nz,1);zeros(Nz,1);zeros(Nz,1)];

% Heater 2 is under lump 22-28
%   Since it spans 7 lumps, the input value is divided by 7
% B_heater2 = [zeros(Nz,1);zeros(Nz,1);zeros(Nz,1);ones(Nz,1)./Nz;zeros(Nz,1);zeros(Nz,1);zeros(Nz,1);zeros(Nz,1);zeros(Nz,1);0];
B_heater2 = [zeros(Nz,1);zeros(Nz,1);zeros(Nz,1);ones(Nz,1)./Nz;zeros(Nz,1);zeros(Nz,1);zeros(Nz,1);zeros(Nz,1);zeros(Nz,1)];

% Heater 3 is under lump 29-34
%   Since it spans 7 lumps, the input value is divided by 7
% B_heater3 = [zeros(Nz,1);zeros(Nz,1);zeros(Nz,1);zeros(Nz,1);ones(Nz,1)./Nz;zeros(Nz,1);zeros(Nz,1);zeros(Nz,1);zeros(Nz,1);0];
B_heater3 = [zeros(Nz,1);zeros(Nz,1);zeros(Nz,1);zeros(Nz,1);ones(Nz,1)./Nz;zeros(Nz,1);zeros(Nz,1);zeros(Nz,1);zeros(Nz,1)];

% Heater 4 is under lump 35-41
%   Since it spans 7 lumps, the input value is divided by 7
% B_heater4 = [zeros(Nz,1);zeros(Nz,1);zeros(Nz,1);zeros(Nz,1);zeros(Nz,1);ones(Nz,1)./Nz;zeros(Nz,1);zeros(Nz,1);zeros(Nz,1);0];
B_heater4 = [zeros(Nz,1);zeros(Nz,1);zeros(Nz,1);zeros(Nz,1);zeros(Nz,1);ones(Nz,1)./Nz;zeros(Nz,1);zeros(Nz,1);zeros(Nz,1)];

% Heater 5 is under lump 42-56
%   Since it spans 14 lumps, the input value is divided by 14
% B_heater5 = [zeros(Nz,1);zeros(Nz,1);zeros(Nz,1);zeros(Nz,1);zeros(Nz,1);zeros(Nz,1);ones(Nz,1)./(2*Nz);ones(Nz,1)./(2*Nz);zeros(Nz,1);0];
B_heater5 = [zeros(Nz,1);zeros(Nz,1);zeros(Nz,1);zeros(Nz,1);zeros(Nz,1);zeros(Nz,1);ones(Nz,1)./(2*Nz);ones(Nz,1)./(2*Nz);zeros(Nz,1)];

% Combine the seperate heater matrices to a single one, with each heater
% its own input --> Resulting in 5 inputs
B_heaters = [zeros((Ny-1)*Nlumps,1),zeros((Ny-1)*Nlumps,1),zeros((Ny-1)*Nlumps,1),zeros((Ny-1)*Nlumps,1),zeros((Ny-1)*Nlumps,1);...
             B_heater1,             B_heater2,             B_heater3,             B_heater4,             B_heater5;...
             zeros(29,1),           zeros(29,1)            zeros(29,1)            zeros(29,1)            zeros(29,1)];

B_heaters = [B_heaters [LTM.GA zeros(1,(Ny-2)*Nlumps+29)]' [zeros(1,(Ny)*Nlumps) waterVelocity*c.mech.density.water*c.therm.specHeat.water zeros(1,28)]'];







