# Thermophysical model for rough surfaces
## Introduction and motivation
The heat flux incident upon the surface of airless planetary bodies is dominated by solar radiation during the day, and by thermal emission from topography at night. The model presented here calculates the temeprature distribution of airless surfaces accounting for insolation, reflected and emitted radiation and subsurface conduction. The model with written in Matlab to emphasize readability and ease of use. 

## Methods
The basic equations used to write this model can here found [here](https://www.sciencedirect.com/science/article/abs/pii/S001910351630656X). 

## Basic instructions and validations
### Installation
Download or clone this repository into a directory of your choosing. The main model directory is composed of three directories and a script file called ```thermodel.m```. The directory ```src``` contains all the source files, ```config``` the config files and ```input``` the topography input file of size NxN (a .mat file with a custom DTM).

After running the main script of the application, two more directories will be created: ```output``` and ```logs```, containing (wait for it...) the output files and logs. ```output``` will contain the main output file, ```Tsurf.mat``` - a NxNxM matlab matrix showing the surface temperature of the topography for time steps 1...M.

### Illumination model
One way to validate this illumination model is to compare its output to that of an analytic model. For example, according to [this](https://www.sciencedirect.com/science/article/abs/pii/001910359290016Z) model, the temperature of a lunar hemispherical (bowl shaped) crater with depth/diameter = 0.2, found at latitude 80 degrees assuming zero obliquity, is ~154 K. To test our code compared to this analytic result, use the following settings:
In file ```config/makeSettings.m```, set the following parameters:
```
    % Simulation settings:
    settings.initializationLengthInSolarDays                     = 0;
    settings.timeStepsPerDayIni                                 = 20;
    settings.simulationLengthInSolarDays                        = 1;
    settings.timeStepsPerDaySim                                 = 6;
    settings.runRCM                                             = true;
    settings.runVFM                                             = true;
    settings.runShadow                                          = true;
    settings.runFluxAndTemperature                              = true;
    settings.calcEquiTemperature                                = true;
    settings.removePreviousFiles                                = false;

    % More settings:
    settings.stopIfErrorsOccur                                  = true;
    settings.createRandomSurface                                = false;
    settings.createSphericalCrater                              = true;
    settings.numberOfScatteringIterations       = 3;
    settings.saveSubsurfTempEveryXSteps         = 3;
    settings.updateTempDepVarXTimesADay         = 20;
    settings.finiteSunArea                      = 9;    % Currently only works for 9 and 1.
    settings.automaticThermalParameters         = true;

    save([simDir,'/config/bin/settings.mat'], 'settings');
```

To change the simulation parameters, edit file ```config/makeConfig.m```:
```
    % Properties related to model-generated surfaces:
    mapProp.mapSize                             = 101;                   % Pixels.
    % Spherical craters:
    mapProp.depthToDiameter                     = 0.2;
    mapProp.sphericalCraterRadius                = 48;                   % Pixels.
    % Random surfaces:
    mapProp.rmsSlope                            = 20;                   % Degrees.

    %% Physical properties:
    % General parameters:
    physProp.distanceFromSun                     = 1;                 % AU.
    physProp.sunRadius                           = 695e3;             % km
    physProp.solarDeclination                    = 0;                 % Degrees.
    physProp.rotationPeriod                      = 86400 * 27.3;       % Seconds.
    physProp.initialSubSurfTemp                  = 150;
    physProp.geothermalFlux                      = 0;                    % W m^2;

    % Optical parameters:
    physProp.albedo                             = 0.136;
    physProp.meanEmissivity                     = 0.95;

```

Note we set ```settings.calcEquiTemperature = true```, and make sure ```mapProp.mapSize``` > 2 * ```mapProp.sphericalCraterRadius```.

Next, make sure matlab is properly [aliased](http://manpages.ubuntu.com/manpages/trusty/man1/alias.1posix.html) in your ```.bashrc``` or ```.bash_profile``` file, and run:

``` 
matlab -nodesktop -nosplash -r thermodel 
```

## Known issues and bugs:
1. The conduction model becomes unstable for temperature dependent heat conductivity and volumetric heat capacity if those are updated at too-large time steps. It is recommended to update these parameters *every* time step when solving problems with highly variable irradiation.
2. ```settings.finiteSunArea``` only works for 1 or 9 and ```settings.stopIfErrorsOccur``` and ```settings.removePreviousFiles``` flags do not work at this stage.
