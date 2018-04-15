%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thermal model of airless planetary bodies      %
% Written by Lior Rubanenko                      %
% Weizmann Institute of Science, Rehovot, Israel %
% University of Los-Angeles, California, USA 	 %
% October 2014 -June 2015                        %
% ...                                            %
% Modified Sep 2016 to support variable thermal  %
% parameters.                                    %
% Modified Jan 2017 to support a finite Sun      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%
% PREPERATIONS %
%%%%%%%%%%%%%%%%
% Clear memory.
clear;

% Store the simulation dir in a variable for later use.
simDir = pwd;

% Add the source directory to the script path:
if exist([simDir,'/src'],'dir')
    addpath([simDir,'/src']);
    addpath([simDir,'/config']);
else
    writeToLog('Cannot continue, src directory is missing.', true);
end


% Read configuation and settings file:
% Remove old binaries
if exist([simDir,'/config/bin'],'dir')
    rmdir([simDir,'/config/bin'],'s');
end

% Read the settings file:
if exist([simDir,'/config/makeSettings.m'],'file')
    % Create the bin in the config directory to store settings in:
    mkdir([simDir, '/config/bin']);
    global settings;
    settings = makeSettings(simDir);
else
    error('Cannot continue, makeConfig.m or makeSettings.m are missing.');
end

% Read the configuration files:
if exist([simDir,'/config/makeConfig.m'],'file')
    [mapProp, physProp, constants] = makeConfig(simDir);
else
    writeToLog('Cannot continue, config/makeConfig.m is missing.', true);
end

% Remove previous log file:
if exist([simDir, '/logs/log.txt'], 'file')
    delete([simDir,'/logs/log.txt']);
end

% Clear previous output files if requested:
if (settings.removePreviousFiles)
    if exist([simDir, '/output'], 'dir')
        rmdir([simDir, '/config/bin'],'s');
        rmdir([simDir, '/output'],'s');
    end
end

% Create output, config, input and logs directories:
directoriesToCreate = {'/input', '/config/bin', '/logs/RCM', '/logs/VFM', '/logs/Flux', '/logs/Shadow', '/logs/Temperature', '/output/Shadow', '/output/buffTotalFlux', '/output/buffTss', '/output/RCM', '/output/VFM', '/output/Shadow', '/output/Flux', '/output/temperature'};

for ii=1:length(directoriesToCreate)
    if ~exist([simDir, directoriesToCreate{ii}], 'dir')
        mkdir([simDir, directoriesToCreate{ii}]);
    end
end

writeToLog('THERMODEL initialized.', true);
writeToLog('Reading configuration files...', true);
writeToLog('Validating input.', true);
validateInput(simDir);
% Update properties after validation:
load([settings.dirPath.config, 'bin/physProp.mat']);
load([settings.dirPath.config, 'bin/mapProp.mat']);

% Write settings to log:
writeToLog('Simulation settings:');
writeToLog(['Initialization length: ',num2str(settings.initializationLengthInSolarDays), ' solar days.'], true);
writeToLog(['Time steps per day (initialization): ',num2str(settings.timeStepsPerDayIni), '.'], true);
writeToLog(['Simulation length: ',num2str(settings.simulationLengthInSolarDays), ' solar days.'], true);
writeToLog(['Time steps per day (Simulation): ',num2str(settings.timeStepsPerDaySim), '.'], true);
writeToLog(['Make ray casting matrix: ',logToStr(settings.runRCM), '.'], true);
writeToLog(['Make view factor matrix: ',logToStr(settings.runVFM), '.'], true);
writeToLog(['Calculate shadows: ',logToStr(settings.runShadow), '.'], true);
writeToLog(['Calculate flux and temperature: ',logToStr(settings.runFluxAndTemperature), '.'], true);
writeToLog(['Calculate only equilibrium temperatures: ',logToStr(settings.calcEquiTemperature), '.'], true);
writeToLog(['Simulation Latitude: ',num2str(mapProp.latitude), '.'], true);


%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERAL CONFIGURATION %
%%%%%%%%%%%%%%%%%%%%%%%%%
% Load the topography (Z) file and handle it.
% If the user has selected a readymade topography:
if settings.createSphericalCrater
    % First, remove the old topography in the input directory:
    if exist([settings.dirPath.input,'Z.mat'], 'file');
        delete([settings.dirPath.input,'Z.mat']);
    end
    writeToLog(['Creating topography: spherical crater with depth to diameter ', num2str(mapProp.depthToDiameter)],true);
    Z = sphericalCrater(mapProp.sphericalCraterRadius, mapProp.depthToDiameter, mapProp.mapSize);
    
elseif settings.createRandomSurface
    % First, remove the old topography in the input directory:
    if exist([settings.dirPath.input,'Z.mat'], 'file');
        delete([settings.dirPath.input,'Z.mat']);
    end
    
    writeToLog(['Creating random topography: Gaussian slope distribution with RMS slope ', num2str(mapProp.rmsSlope)], true);
    Z = gaussianRandomSurface(mapProp.mapSize, mapProp.rmsSlope);
else
    writeToLog('Using a custom topography from the input folder.', true);
    if exist([settings.dirPath.input,'Z.mat'],'file')
        load([settings.dirPath.input,'Z.mat']);
    else
        writeToLog('Cannot continue, Z.mat (topography) is missing from input folder.', true);
    end
end

% Convert to natural units ("pixels") by diving by the scaling factor:
Z = Z / mapProp.scaleFactor; save([settings.dirPath.input,'Z.mat'], 'Z');

% Get the map size from Z.
mapProp.mapSize = length(Z); save([simDir,'/config/bin/mapProp.mat'], 'mapProp');

% Read the data file (azimuth, incidence and time):
if exist([simDir,'/config/data.m'],'file')
    run([simDir,'/config/data.m']);
    
    % Calculate the time vector in case no time vector has been given:
    if isempty(time)
        time = calcTimeVector(solarAzimuth, solarZenithAngle, simDir);
    end
else
    writeToLog('A custom config/data.m file was not found. Creating data.m from input parameters.', true);
    
    if ( settings.initializationLengthInSolarDays > 0 && settings.timeStepsPerDayIni > 0)
        [solarAzimuthIni, solarZenithAngleIni, timeIni] = calcSunTrajectory(physProp.solarDeclination, mapProp.latitude, settings.timeStepsPerDayIni, settings.initializationLengthInSolarDays, physProp.rotationPeriod);
    else
        solarAzimuthIni = [];
        solarZenithAngleIni = [];
        timeIni = 0;
    end
    
    [solarAzimuthSim, solarZenithAngleSim, timeSim] = calcSunTrajectory(physProp.solarDeclination, mapProp.latitude, settings.timeStepsPerDaySim, settings.simulationLengthInSolarDays ,physProp.rotationPeriod);
    
    % Combining to one vector:
    solarAzimuth = [solarAzimuthIni solarAzimuthSim];
    solarZenithAngle = [solarZenithAngleIni solarZenithAngleSim];
    % In order to get the correct time vector, one must calculate timeSim+the
    % last time step in timeIni plus the new dt.
    time = [timeIni timeSim+timeIni(end)+(timeSim(2)-timeSim(1))];
    
    % Since the time vector has an additional value at the beginning (zero) and due to the CN method requiring one extra time step, the solar azimuth and incidence vectors
    % will be extended:
    solarZenithAngle = [solarZenithAngle solarZenithAngleSim(1:(length(time)-length(solarZenithAngle)))];
    solarAzimuth = [solarAzimuth solarAzimuthSim(1:(length(time)-length(solarAzimuth)))];
    
end

writeToLog(['Simulation solar Azimuth Vector:', toCommaDelimitedString(solarAzimuth)]);
writeToLog(['Simulation solar zenith angle Vector:', toCommaDelimitedString(solarZenithAngle)]);
writeToLog(['Simulation solar time:', toCommaDelimitedString(time)]);


% Save to the config/bin directory:
save([simDir,'/config/bin/solarAzimuth.mat'], 'solarAzimuth');
save([simDir,'/config/bin/solarZenithAngle.mat'], 'solarZenithAngle');
save([simDir,'/config/bin/time.mat'], 'time');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RAY CASTING MATRIX (RCM) CALCULATION %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (settings.runRCM)
    if (exist('./output/RCM/RCM1.mat','file'))
        writeToLog('Ray-casting matrix already exists. Skipping calculation.',true);
    else
        writeToLog('Calculating ray-casting matrix.',true);
        createRCM(1:numel(Z), simDir);    
        writeToLog('Finished RCM calculation.');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VIEW FACTOR MATRIX (VFM) CALCULATION %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (settings.runVFM)
    if (exist('./output/VFM/VFM1.mat','file'))
        writeToLog('View factor matrix already exists. Skipping calculation.',true);
    else
        writeToLog('Calculating view factor matrix.',true);
        % Just run the view factor script to all elements in Z.
        createVFM(1:numel(Z), simDir);
        writeToLog('Finished VFM calculation.');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SHADOW AND SOLAR FLUX CALCULATION %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
writeToLog('Begining shadow and solar flux calculation.');

if (settings.runShadow)
    writeToLog('Calculating shadow matrices.',true);
    fprintf('Progress: ');
    for timeStep=1:length(solarZenithAngle)       
        % Show progress:
        if (mod(timeStep,1) == 0)
            progressBar = [num2str(round(timeStep/length(solarZenithAngle) * 100,2)), '%%'];
            fprintf(progressBar);
            fprintf(repmat('\b',1,length(progressBar)-1))
        end
        
        % Quit the script if there is already a solar flux matrix for the current timestep:
        if exist([settings.dirPath.output, 'Shadow/solarFluxMatrix_', num2str(solarAzimuth(timeStep)), '_', num2str(solarZenithAngle(timeStep)) ,'.mat'], 'file')
            writeToLog(['Solar flux map for incidence ', num2str(solarZenithAngle(timeStep)), ' and azimuth ', num2str(solarAzimuth(timeStep)), ' already exists.']);
            continue;
        end
        
        % Calculate the Sun's angular size in the sky:
        [finiteSunZenithAngle, finiteSunAzimuth] = finiteSunCoordinates(solarZenithAngle(timeStep), solarAzimuth(timeStep), simDir);
        
        % If the solar incidence matrix is not in < 90 degrees + angular solar radius, the sun is below the horizon (everything is shadowed):
        for sunPix = 1:numel(finiteSunZenithAngle)
            if (finiteSunZenithAngle(sunPix) >= 90)
                shadowMatrix(:,:,sunPix) = ones(size(Z));
                shadowDepthMatrix(:,:,sunPix) = NaN(size(Z));
                solarFluxMatrix(:,:,sunPix) = zeros(size(Z));
            else
                [buffShad, buffShadDepth, buffSolar] = calcShadowMatrix(finiteSunZenithAngle(sunPix), finiteSunAzimuth(sunPix), simDir);
                shadowMatrix(:,:,sunPix) = buffShad;
                shadowDepthMatrix(:,:,sunPix) = buffShadDepth;
                solarFluxMatrix(:,:,sunPix) = buffSolar;
            end
        end
        
        shadowMatrix = all(shadowMatrix,3);
        shadowDepthMatrix = min(shadowDepthMatrix,[],3);
        solarFluxMatrix = sum(solarFluxMatrix,3);
        
        save([settings.dirPath.output, 'Shadow/shadowMatrix_', num2str(solarAzimuth(timeStep)), '_', num2str(solarZenithAngle(timeStep)) ,'.mat'], 'shadowMatrix');
        save([settings.dirPath.output, 'Shadow/shadowDepthMatrix_', num2str(solarAzimuth(timeStep)), '_', num2str(solarZenithAngle(timeStep)) ,'.mat'], 'shadowDepthMatrix');
        save([settings.dirPath.output, 'Shadow/solarFluxMatrix_', num2str(solarAzimuth(timeStep)), '_', num2str(solarZenithAngle(timeStep)) ,'.mat'], 'solarFluxMatrix');
    end
    
    % Calculate permanent shadow:
    if ~exist([settings.dirPath.output, 'Shadow/permShadowMatrix.mat'], 'file')
        tt = 1;
        
        % Create permanent shadow (and depth) matrices
        for timeStep=1:length(solarZenithAngle)
            if (solarZenithAngle(timeStep) < 90)
                % Load the shadow and shadow depth matrices:
                load([settings.dirPath.output, 'Shadow/shadowMatrix_', num2str(solarAzimuth(timeStep)), '_', num2str(solarZenithAngle(timeStep)),'.mat']);
                load([settings.dirPath.output, 'Shadow/shadowDepthMatrix_', num2str(solarAzimuth(timeStep)), '_', num2str(solarZenithAngle(timeStep)),'.mat']);
                
                % Save the shadow matrices to a 3D array in order to form the
                % permanent shadow (and depth) matrices:
                buffPermShadowMatrix(:,:,tt) = shadowMatrix;
                buffPermShadowDepthMatrix(:,:,tt) = shadowDepthMatrix;
                tt = tt + 1;
            end
        end
        
        permShadowMatrix = all(buffPermShadowMatrix, 3);
        permShadowDepthMatrix = min(buffPermShadowDepthMatrix, [], 3);
        
        save([settings.dirPath.output, 'Shadow/permShadowMatrix.mat'], 'permShadowMatrix');
        save([settings.dirPath.output, 'Shadow/permShadowDepthMatrix.mat'], 'permShadowDepthMatrix');
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TOTAL FLUX AND TEMPERATURE CALCULATION %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
writeToLog('Starting total flux and temperature calculation.');

% Create all required data output files:
writeToLog('Preallocating memory to all output variables.');
if (~exist([settings.dirPath.output,'Tsurf.mat'],'file'))
    Tsurf=zeros(length(Z), length(Z), length(time));
    % Save Tsurf as a -v7.3 file so that it could be partially loaded into memory:
    save([settings.dirPath.output,'Tsurf.mat'],'Tsurf','-v7.3');

    if (exist([settings.dirPath.output,'Tsurf.mat'], 'file'))
        writeToLog('Tsurf.mat successfully preallocated.');
    else
        writeToLog('Tsurf.mat was not successfully preallocated.');
    end

    totalFlux = zeros(size(Tsurf));
    % Save Tsurf as a -v7.3 file so that it could be partially loaded into memory:
    save([settings.dirPath.output,'totalFlux.mat'],'totalFlux','-v7.3');

    if (exist([settings.dirPath.output,'totalFlux.mat'], 'file'))
        writeToLog('totalFlux.mat successfully preallocated.');
    else
        writeToLog('totalFlux.mat was not successfully preallocated.');
    end
    clear totalFlux;

    emissionFlux = zeros(size(Tsurf));
    % Save Tsurf as a -v7.3 file so that it could be partially loaded into memory:
    save([settings.dirPath.output,'emissionFlux.mat'],'emissionFlux','-v7.3');

    if (exist([settings.dirPath.output,'emissionFlux.mat'], 'file'))
        writeToLog('emissionFlux.mat successfully preallocated.');
    else
        writeToLog('emissionFlux.mat was not successfully preallocated.');
    end

    clear emissionFlux;

    reflectionFlux = zeros(size(Tsurf));
    % Save Tsurf as a -v7.3 file so that it could be partially loaded into memory:
    save([settings.dirPath.output,'reflectionFlux.mat'],'reflectionFlux','-v7.3');

    if (exist([settings.dirPath.output,'reflectionFlux.mat'], 'file'))
        writeToLog('reflectionFlux.mat successfully preallocated.');
    else
        writeToLog('reflectionFlux.mat was not successfully preallocated.');
    end

    clear reflectionFlux;

else
    writeToLog('Old Tsurf.mat was found in output folder, continuing previous simulation.', 'true');
    load([settings.dirPath.output,'Tsurf.mat']);
end

% Find the first zero element in the kth dimension of Tsurf to use as the first timestep.
firstk = find(squeeze(all(all(~Tsurf(:,:,2:end)))),1);

[~,~,lastk] = size(Tsurf);
writeToLog(['Simulation will begin from time step ', num2str(firstk), ' and end at time step ', num2str(lastk)]);
% Clear Tsurf to save memory:
clear Tsurf;

% If flux and temperature calculation was requested:
if (settings.runFluxAndTemperature)
    writeToLog('Beginning flux and temperature calculation');

    % If the model runs as a stand alone:
    for (timeStep = firstk:(lastk-1))
        if (mod(timeStep, 5)) == 0
            copyfile('output/Tsurf.mat','output/Tsurf_bu.mat', 'f');
        end
        
        writeToLog(['Beginning time step ', num2str(timeStep), '.'], true);
        
        for scatteringIteration = 1:settings.numberOfScatteringIterations
            writeToLog(['Beginning scattering event ', num2str(scatteringIteration), '.'], true);
            calcIncomingFlux(1:numel(Z), timeStep, scatteringIteration, simDir);
            jobSize = mapProp.mapSize.^2; lastJobSize = 0;
            
            % This is temporary, until 
            numberOfJobs = 1;
            consolidateFlux(numberOfJobs, jobSize, lastJobSize, timeStep, simDir);
        end
        
        if (settings.calcEquiTemperature == false)
            writeToLog('Calculating subsurface conduction.',true);
            calcTemperature(1:numel(Z), time(timeStep + 1) - time(timeStep), timeStep, simDir);
            consolidateTemperature(settings.numberOfJobs, jobSize, lastJobSize, timeStep, simDir);
        else
            writeToLog('Calculating only equilibrium temperatures',true);
        end
        
        writeToLog(['TimeStep ', num2str(timeStep), ' has been completed.'], true);
    end
end

writeToLog('Simulation has ended.',true);
