function[F] = consolidateFlux(numberOfJobs, jobSize, lastJobSize, timeStep, simDir)
writeToLog('Consolidating total flux into totalFlux.mat.');

% Load configuration files and assign variables:
load([simDir, '/config/bin/settings.mat']);
load([settings.dirPath.config, 'bin/constants.mat']);
load([settings.dirPath.config, 'bin/physProp.mat']);
load([settings.dirPath.config, 'bin/mapProp.mat']);

fTsurf = matfile([settings.dirPath.output,'Tsurf.mat'],'writable',true);
fTotFlux = matfile([settings.dirPath.output,'totalFlux.mat'],'writable',true);
fEmis = matfile([settings.dirPath.output,'emissionFlux.mat'],'writable',true);
fRef = matfile([settings.dirPath.output,'reflectionFlux.mat'],'writable',true);

totalFlux = zeros(mapProp.mapSize);
emissionFlux = zeros(mapProp.mapSize);
reflectionFlux = zeros(mapProp.mapSize);

% Take the first cell of each subsurface temperature vector and put it in totalFlux:
for ii=1:numberOfJobs
	if (numberOfJobs > 1)
		[elementRangeMax, elementRangeMin] = calculateIterationSize(ii, jobSize, lastJobSize, numberOfJobs, mapProp.mapSize);
       	elementRange = elementRangeMin:elementRangeMax;

    elseif (numberOfJobs == 1)
    	elementRange = 1:mapProp.mapSize^2;
    end

    % Load the correct flux range in order to consolidate it:
	load([settings.dirPath.output, 'buffTotalFlux/buffTotalFlux_', num2str(elementRange(1)),'.mat']);
	load([settings.dirPath.output, 'buffTotalFlux/bufEmission_', num2str(elementRange(1)),'.mat']);
    load([settings.dirPath.output, 'buffTotalFlux/bufVisReflection_', num2str(elementRange(1)),'.mat']);
    

	for jj=1:length(elementRange)
        totalFlux(elementRange(jj)) = buffTotalFlux(jj);
		emissionFlux(elementRange(jj)) = bufEmission(jj);
		reflectionFlux(elementRange(jj)) = bufVisReflection(jj);
	end

end

fTotFlux.totalFlux(:,:,timeStep) = totalFlux;
fEmis.emissionFlux(:,:,timeStep) = emissionFlux;
fRef.reflectionFlux(:,:,timeStep) = reflectionFlux;

% Calculate the equilibrium (or "quasi equilibrium" in case conduction is
% on) temperature:
fTsurf.Tsurf(:,:,timeStep) = nthroot(totalFlux/(constants.sb .* physProp.meanEmissivity), 4);

