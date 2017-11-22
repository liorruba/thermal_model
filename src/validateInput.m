function[] = validateInput(simDir)
% Load configuration files and assign variables:
load([simDir, '/config/bin/settings.mat']);
load([settings.dirPath.config, 'bin/physProp.mat']);
load([settings.dirPath.config, 'bin/mapProp.mat']);


if mapProp.latitude > 90 || mapProp.latitude < 0
    writeToLog('ERROR: latitude (mapProp.latitude) must be between 0 and 90');
    error('northmost latitude (mapProp.latitude) must be between 0 and 90');
end

if mapProp.zAxis(1) ~= 0
    writeToLog('The first element of the z axis (mapProp.zAxis) must be 0.');
    error('The first element of the z axis (mapProp.zAxis) must be 0.');
end

% Force mapProp.zAxis to be a cloumn vector:
if isrow(mapProp.zAxis)
    mapProp.zAxis = mapProp.zAxis'; 
end

if mapProp.mapSize < 0
    writeToLog('Map size (mapProp.mapSize) must be > 0');
    error('Map size (mapProp.mapSize) must be > 0');
end

if (mapProp.depthToDiameter > 1 || mapProp.depthToDiameter < 0) && settings.createSphericalCrater	== true
    writeToLog('Depth/diameter ratio must be between 0 and 1.');
    error('Depth/diameter ratio must be between 0 and 1.');
end 

if mapProp.sphericalCraterRadius >= mapProp.mapSize/2 && settings.createSphericalCrater	== true
    writeToLog('The radius of the spherical crater must be > 1/2 map size');
    error('The radius of the spherical crater must be > 1/2 map size');
end 
    
if settings.automaticThermalParameters == false
    if isscalar(physProp.specificHeatCapacity)
            physProp.specificHeatCapacity = ones(length(mapProp.zAxis)-1,1).* physProp.specificHeatCapacity;
    end

    if isscalar(physProp.density)
        physProp.density = ones(length(mapProp.zAxis)-1,1) .* physProp.density;
    end

    physProp.thermalConductivity = physProp.thermalInertia.^2 ./ physProp.specificHeatCapacity ./ physProp.density;
    
end

if (settings.finiteSunArea ~= 9) && (settings.finiteSunArea ~= 1)
   error('settings.finiteSunArea must be either 1 (point Sun) or 9 (finite Sun).');
end

save([settings.dirPath.config,'bin/physProp.mat'], 'physProp');
save([settings.dirPath.config,'bin/mapProp.mat'], 'mapProp');
end