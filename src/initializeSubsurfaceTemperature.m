function [Tsurf] = initializeSubsurfaceTemperature(simDir)
% Load configuration files and assign variables:
load([simDir, '/config/bin/settings.mat']);
load([settings.dirPath.config, 'bin/constants.mat']);
load([settings.dirPath.config, 'bin/physProp.mat']);
load([settings.dirPath.config, 'bin/mapProp.mat']);

writeToLog('Initializing surface temperatures.')

dirShad = [settings.dirPath.output,'Shadow/'];

fileList = dir([dirShad,'solarFlux*']);

ii=1;
for file=fileList'
    load([dirShad, file.name]);
    
    s(:,:,ii) = nthroot((1 - physProp.albedo) * solarFluxMatrix / (physProp.meanEmissivity * constants.sb),4);
    ii = ii + 1;
end

Tsurf = sum(s,3)./ii;
Tsurf(Tsurf<100) = 100;
