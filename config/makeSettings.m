function [settings] = makeSettings(simDir)
    %
    % Main settings file:
    %

    % File paths:
    settings.dirPath.src = [simDir,'/src/'];
    settings.dirPath.output = [simDir,'/output/'];
    settings.dirPath.input = [simDir,'/input/'];
    settings.dirPath.logs = [simDir,'/logs/'];
    settings.dirPath.config = [simDir,'/config/'];

    % In order to parallelize the model, enter how many jobs to run on each iteration:
    % (FOR THE STANDALONE VERSION, SET THIS TO 1.)
    settings.numberOfJobs 					= 1;

    % Simulation name:
    [~,simDirName,~] = fileparts(simDir);
    settings.simulationName 					= simDirName;

    % Simulation settings:
    settings.initializationLengthInSolarDays	= 0;
    settings.timeStepsPerDayIni					= 20;
    settings.simulationLengthInSolarDays		= 1;
    settings.timeStepsPerDaySim					= 720;
    settings.runRCM								= false;
    settings.runVFM								= false;
    settings.runShadow							= true;
    settings.runFluxAndTemperature				= false;
    settings.calcEquiTemperature				= false;
    settings.removePreviousFiles				= false;

    % More settings:
    settings.stopIfErrorsOccur					= true;
    settings.createRandomSurface				= true;
    settings.createSphericalCrater				= false;
    settings.numberOfScatteringIterations       = 3;
    settings.saveSubsurfTempEveryXSteps         = 3;
    settings.updateTempDepVarXTimesADay         = 20;
    settings.finiteSunArea                      = 9;    % Currently only works for 9 and 1.
    settings.automaticThermalParameters         = true;

    save([simDir,'/config/bin/settings.mat'], 'settings');
end
