function[] = calcTemperature(elementRange, dt, timeStep, simDir)
% Load configuration files and assign variables:
load([simDir, '/config/bin/settings.mat']);
load([settings.dirPath.config, 'bin/constants.mat']);
load([settings.dirPath.config, 'bin/physProp.mat']);
load([settings.dirPath.config, 'bin/mapProp.mat']);

writeToLog('Calculating subsurface temperatures.');

z = mapProp.zAxis;

% Load the incoming flux matrix:
fTsurf = matfile([settings.dirPath.output,'Tsurf.mat'],'writable',true);
fTotalFlux = matfile([settings.dirPath.output,'totalFlux.mat'],'writable',false);

% Load the correct subsurface temperature matrix.
% Only during the equilibration process, set the subsurface temperature to the mean diurnal surface temperature of the previous day:
if ( mod(timeStep, settings.timeStepsPerDayIni) == 1 && timeStep ~= 1 && timeStep < (settings.timeStepsPerDayIni * settings.initializationLengthInSolarDays))
    TMean = mean(fTsurf.Tsurf(:,:,timeStep-settings.timeStepsPerDayIni:timeStep-1), 3);
    T0 = bsxfun(@times, ones(length(z), length(elementRange)), TMean(elementRange));
else
    % During the first time step, set the subsurface temperature to the equilibrium temperature of the surface:
    if ( timeStep == 1 )
        Tsurf0 = physProp.initialSubSurfTemp * ones(mapProp.mapSize);
        T0 = bsxfun(@times, ones(length(z), length(elementRange)),Tsurf0(elementRange));
    else
        % During any other time step, load the last step subsurface temperature:
        load([settings.dirPath.output, 'buffTss/T0_', num2str(elementRange(1)),'.mat']);
    end
end

% Iterate all facets for the current time step:
for facetLinearIndex = 1:length(elementRange)
    totalFlux = fTotalFlux.totalFlux(:,:,timeStep);
    
    if (timeStep > 1)
        lastStepTotalFlux = fTotalFlux.totalFlux(:,:,timeStep-1);
    else
        lastStepTotalFlux = zeros(size(fTotalFlux.totalFlux(:,:,timeStep)));
    end
    
    incomingRadiation = totalFlux(elementRange(facetLinearIndex));
    lastStepIncomingRadiation = lastStepTotalFlux(elementRange(facetLinearIndex));
    
    % If the specific heat capacity is not provided by the user,
    % update it during the day (change via the settings file), according to
    % the appendix of Hayne et al., unpublished.
    if settings.automaticThermalParameters == true
        if mod(timeStep, ceil(settings.timeStepsPerDayIni/settings.updateTempDepVarXTimesADay)) == 1 || mod(timeStep, ceil(settings.timeStepsPerDaySim/settings.updateTempDepVarXTimesADay)) == 1 ||...
                ~isfield(physProp,'thermalConductivity') || ~isfield(physProp,'density') || ~isfield(physProp,'specificHeatCapacity')
                
            % Pass just the subsurface (T0(2:end)) to update the thermal
            % parameters:
            % Log every new iteration:
            if facetLinearIndex == 1
                writeToLog('Updating thermal parameters.',true);
            end
            
            updateThermalParameters(settings, T0(2:end,facetLinearIndex));
            % Reload the updated properties file:
            load([settings.dirPath.config, 'bin/physProp.mat']);
        end
    end
    
    T1(:, facetLinearIndex) = heatDiffusionSemiImplicit_1d(z, dt, T0(:,facetLinearIndex), physProp.specificHeatCapacity, physProp.density, physProp.thermalConductivity, physProp.meanEmissivity, incomingRadiation, physProp.geothermalFlux);
    
    % Fix any temperature instability that is due to an abrupt change in
    % the incoming flux (due to the Crank-Nicolson method):
    if T1(1, facetLinearIndex) - T0(1, facetLinearIndex) > 30 && abs(lastStepIncomingRadiation - incomingRadiation) > 10
        
        T1(:, facetLinearIndex) = fixTemperatureInstability(z, dt, T0(:,facetLinearIndex), physProp.specificHeatCapacity, physProp.density, physProp.thermalConductivity, physProp.meanEmissivity, lastStepIncomingRadiation, incomingRadiation, physProp.geothermalFlux,2);
    end
    
end

% Changing the name back to T0 so that the .mat files loads with the correct variable name:
T0 = T1;

% Every x time steps save the subsurface temperatures to a new file:
if ( mod(timeStep, settings.saveSubsurfTempEveryXSteps) == 0 )
    save([settings.dirPath.output, 'TsubSurf',num2str(timeStep),'_',num2str(elementRange(1)),'.mat'], 'T0');
end

% Saving to a temporary file to be consolidated later:
save([settings.dirPath.output, 'buffTss/T0_', num2str(elementRange(1)),'.mat'], 'T0');

% Open the lock file and append the char "1":
incrementLock(simDir);
