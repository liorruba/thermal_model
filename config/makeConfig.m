function [mapProp, physProp, constants] = makeConfig(simDir)
    %
    % Physical and geographical settings:
    %

    %% Map properties:
    mapProp.scaleFactor                         = 1;                    % m pixels^-1
    mapProp.latitude                    = 80;                   % Degrees.
    mapProp.zAxis = 0.01*(cumprod(ones(1,36) * 1.13) - 1.13);          % Meters

    % Properties related to model-generated surfaces:
    mapProp.mapSize                             = 25;                   % Pixels.
    % Spherical craters:
    mapProp.depthToDiameter                     = 0.2;
    mapProp.sphericalCraterRadius       = 48;                   % Pixels.
    % Random surfaces:
    mapProp.rmsSlope                            = 20;                   % Degrees.

    %% Physical properties:
    % General parameters:
    physProp.distanceFromSun            = 1;                 % AU.
    physProp.sunRadius              = 695e3;        % km
    physProp.solarDeclination           = 0;                  % Degrees.
    physProp.rotationPeriod                     = 86400 * 27.3;      % Seconds.
    physProp.initialSubSurfTemp         = 150;
    physProp.geothermalFlux                     = 0;                    % W m^2;

    % Optical parameters:
    physProp.albedo                             = 0.136;
    physProp.meanEmissivity             = 0.95;

    % Thermal parameters:
    % If settings.autoThermalParameters=true, must provide additional parameters:
    physProp.radiativeConductivityParameter     = 2.7;
    physProp.radiativeConductivityFactor        = physProp.radiativeConductivityParameter / 350^3;
    physProp.surfaceContactConductivity         = 3e-4;             % W m^-1 K^-1.
    physProp.deepContactConductivity            = 5e-3;             % W m^-1 K^-1.
    physProp.Hparameter                         = 0.06;             % m
    physProp.surfaceDensity                     = 1100;             % kg m^-3
    physProp.depthDensity                       = 1800;             % kg m^-3

    % If settings.autoThermalParameters=false, must provide the thermal
    % parameters manually:
    physProp.specificHeatCapacity   = 300;              % J kg^-1 K^-1
    physProp.thermalInertia                     = 25;               % J m^-2 K^-1 s^-0.5
    physProp.density                = 1000;             % J kg^-1 K^-1

    % Constants:
    constants.sb                                        = 5.6704e-8;        % SI
    constants.earthSolarConstant        = 1367;                 % W m^-2

    save([simDir,'/config/bin/constants.mat'], 'constants');
    save([simDir,'/config/bin/physProp.mat'], 'physProp');
    save([simDir,'/config/bin/mapProp.mat'], 'mapProp');
end
                                                                                                                                                              
                                                                                                                                                              
                                                                                                                                                              
                                                                                                                                                              
                                                                              
