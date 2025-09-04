function constants = constsGen(varargin)

vers = {'thermal','mech','elec','chem'};

constants.gen.g = 9.81;

if ~isempty(varargin)
    % disp(vers(nonzeros(double(ismember(vers,varargin)).*double(ones(size(vers))*triu(ones(size(vers,2)))))'))
    toLoadConsts = vers(nonzeros(double(ismember(vers,varargin)).*double(ones(size(vers))*triu(ones(size(vers,2)))))');

    if sum(ismember(toLoadConsts,'mech'))
    % Specific heat constants
        % Gotten from - https://www.ecs.csun.edu/~nhuttho/me584/Chapter%206%20Thermal%20Systems.pdf, page 17
        constants.mech.density.alu    = 2702;
        constants.mech.density.copper = 8933;
        constants.mech.density.steel  = 7870;
        constants.mech.density.silver = 10500;
        constants.mech.density.water  = 1000;
    end
    if sum(ismember(toLoadConsts,'thermal'))
    % Specific heat constants
        % Gotten from - https://www.ecs.csun.edu/~nhuttho/me584/Chapter%206%20Thermal%20Systems.pdf, page 13
        % constants.therm.specHeat.alu      = 900/15;
        % constants.therm.specHeat.alu      = 900/2;
        % constants.therm.specHeat.alu      = 900;
        constants.therm.specHeat.alu      = 1300;
        constants.therm.specHeat.copper   = 380;
        constants.therm.specHeat.iron     = 450;
        constants.therm.specHeat.lead     = 130;
        constants.therm.specHeat.silver   = 240;
        constants.therm.specHeat.stone    = 920;
        constants.therm.specHeat.concrete = 960;
        constants.therm.specHeat.glas     = 840;
        constants.therm.specHeat.water    = 4180;
        constants.therm.specHeat.air      = 1000;
        constants.therm.specHeat.RVS      = 500;

    % Heatconduction constants
        % Gotten from - Inleiding Warmte en Stroming, 4B260 & 4B270, page 126
        % constants.therm.heatCond.alu      = 237;
        constants.therm.heatCond.alu      = 180;
        constants.therm.heatCond.copper   = 403;
        constants.therm.heatCond.gold     = 318;
        constants.therm.heatCond.iron     = 80.4;
        constants.therm.heatCond.lead     = 35.5;
        constants.therm.heatCond.mercury  = 8.3;
        constants.therm.heatCond.steel    = 43;
        constants.therm.heatCond.RVS      = 25;
        constants.therm.heatCond.silver   = 429;
        constants.therm.heatCond.stone    = 0.8;
        constants.therm.heatCond.concrete = 0.15;
        constants.therm.heatCond.glas     = 0.78;
        constants.therm.heatCond.wood1    = 0.08;
        constants.therm.heatCond.wood2    = 0.16;
        constants.therm.heatCond.water    = 0.56;
        constants.therm.heatCond.air      = 0.026;
    
    % Convection constants
        % Gotten from [1] - Inleiding Warmte en Stroming, 4B260 & 4B270, page 129
        % Gotten from [2] - https://www.ecs.csun.edu/~nhuttho/me584/Chapter%206%20Thermal%20Systems.pdf, page 11
        constants.therm.Conv.airNatural.min  = 2.5;     % [1]
        constants.therm.Conv.airNatural.nom  = 9.6;      % [1]
        % constants.therm.Conv.airNatural.nom  = 12;      % [1]
        % constants.therm.Conv.airNatural.nom  = 15;      % [1]
        constants.therm.Conv.airNatural.max  = 25;      % [1]
        constants.therm.Conv.airForced.min   = 10;      % [1]
        constants.therm.Conv.airForced.nom   = 200;     % [1]
        constants.therm.Conv.airForced.max   = 500;     % [1]
        constants.therm.Conv.fluidForced.min = 100;     % [1]
        constants.therm.Conv.fluidForced.max = 15000;   % [1]
        constants.therm.Conv.waterFree.min   = 60;      % [2]
        constants.therm.Conv.waterFree.nom   = 5;       % [2]
        constants.therm.Conv.waterFree.max   = 300;     % [2]
        constants.therm.Conv.waterForced.min = 300;     % [2]
        constants.therm.Conv.waterForced.nom = 350;     % [2]
        constants.therm.Conv.waterForced.max = 6000;    % [2]
        constants.therm.Conv.boilWater.min   = 2500;    % [1]
        constants.therm.Conv.boilWater.max   = 25000;   % [1]
    
    % Stefan-Boltzmann constant
        constants.therm.stefBoltz = 5.67e-8; %[Wm^-2K^-4]

    % Heat emmision coefficients
        % Gotten from - Inleiding Warmte en Stroming, 4B260 & 4B270, page 131
        constants.therm.heatEmmision.alu         = 0.09;
        constants.therm.heatEmmision.iron_rust   = 0.74;
        constants.therm.heatEmmision.iron_melted = 0.88;
        constants.therm.heatEmmision.water       = 0.96;
        constants.therm.heatEmmision.skin        = 0.97;
    end

end
end