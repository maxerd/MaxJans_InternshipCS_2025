
% Initialize the E matrix
% E = zeros(Ny*Nlumps+1);
E = zeros(Ny*Nlumps+29);

% Make the E matrix (Determine the thermal masses)
%   Only has diagonal terms with m*c for the respective lump
for i=1:Nlumps
    E(i,i)               = TMA(i).mass/1000*c.therm.specHeat.alu; % For multilayer system
    E(i+Nlumps,i+Nlumps) = TMB(i).mass/1000*c.therm.specHeat.alu; % For multilayer system

    %%%%%%%%%%%%%%%%%% < Semi-random value > %%%%%%%%%%%%%%%%%%
    % Heater thermal mass calculation
    % E(i+2*Nlumps,i+2*Nlumps) = 3.15/1000*c.therm.specHeat.alu; % For multilayer system
    E(i+2*Nlumps,i+2*Nlumps) = 8.5/1000*c.therm.specHeat.alu; % For multilayer system
    % E(i+2*Nlumps,i+2*Nlumps) = 15/1000*c.therm.specHeat.alu; % For multilayer system
    %%%%%%%%%%%%%%%%%% ^ Semi-random value ^ %%%%%%%%%%%%%%%%%%

end

for i=1:29
    % E((end-29+i),(end-29+i)) = 0.0025*c.mech.density.water*c.therm.specHeat.water;
    E((end-29+i),(end-29+i)) = 0.0025*c.therm.specHeat.water;
end

% Define the ambient air thermal mass
% E(end,end) = LAMB.mass/1000*c.therm.specHeat.air;

% Get a vector with the values of the E matrix, used for debugging only
% for i=1:Ny*Nlumps+1
for i=1:Ny*Nlumps
    E_vals(i) = E(i,i);
end