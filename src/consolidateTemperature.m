function[] = consolidateTemperature(numberOfJobs, jobSize, lastJobSize, timeStep, simDir)
% Load configuration files and assign variables:
load([simDir, '/config/bin/settings.mat']);
load([settings.dirPath.config, 'bin/constants.mat']);
load([settings.dirPath.config, 'bin/physProp.mat']);
load([settings.dirPath.config, 'bin/mapProp.mat']);

writeToLog('Consolidating temperatures into Tsurf.mat.');

% Load data files:
load([settings.dirPath.config, 'bin/solarAzimuth.mat']);
load([settings.dirPath.config, 'bin/solarZenithAngle.mat']);

fTsurf = matfile([settings.dirPath.output,'Tsurf.mat'],'writable',true);
Tsurf = zeros(mapProp.mapSize);

% Take the first cell of each subsurface temperature vector and put it in Tsurf:
for ii=1:numberOfJobs
    if (numberOfJobs > 1)
        [elementRangeMax, elementRangeMin] = calculateIterationSize(ii, jobSize, lastJobSize, numberOfJobs, mapProp.mapSize);
        elementRange = elementRangeMin:elementRangeMax;

    elseif (numberOfJobs == 1)
        elementRange = 1:mapProp.mapSize.^2;
    end

    % Load the correct temperature range in order to consolidate it:
    load([settings.dirPath.output, 'buffTss/T0_', num2str(elementRange(1)),'.mat']);   

    for jj=1:length(elementRange)   
        Tsurf(elementRange(jj)) = T0(1,jj);
    end

end

fTsurf.Tsurf(:,:,timeStep) = Tsurf;
