function[] = calcIncomingFlux(elementRange, timeStep, scatteringEvent, simDir)
% Load configuration files and assign variables:
load([simDir, '/config/bin/settings.mat']);
load([settings.dirPath.config, 'bin/constants.mat']);
load([settings.dirPath.config, 'bin/physProp.mat']);
load([settings.dirPath.config, 'bin/mapProp.mat']);

% Load data files:
load([simDir,'/config/bin/solarAzimuth.mat']);
load([simDir,'/config/bin/solarZenithAngle.mat']);

% Load the current step solar matrix:
load([settings.dirPath.output, 'Shadow/solarFluxMatrix_', num2str(solarAzimuth(timeStep)), '_', num2str(solarZenithAngle(timeStep)),'.mat']);

% Initialize a vector in which the temporary values of the solar flux will be saved at. Later, this vector will be consolidated to a full flux matrix:
buffTotalFlux = zeros(size(elementRange));
bufEmission = zeros(size(elementRange));
bufVisReflection = zeros(size(elementRange));

% Load previous step matrices as matfiles:
fRef = matfile([settings.dirPath.output,'reflectionFlux.mat'],'writable',false);
fEmis = matfile([settings.dirPath.output,'emissionFlux.mat'],'writable',true);
fTsurf = matfile([settings.dirPath.output,'Tsurf.mat'],'writable',true);

% Set the first time step to be the equilibrium temperature of the surface:
if scatteringEvent == 1
    TsurfEqui = nthroot(( 1 - physProp.albedo ) .* solarFluxMatrix / (constants.sb .* physProp.meanEmissivity),4);
    
else
    TsurfEqui = fTsurf.Tsurf(:,:,timeStep);
   
    % At night, load the previous time step (to avoid convergence issues):
    if (TsurfEqui == 0)
        TsurfEqui = fTsurf.Tsurf(:,:,timeStep - 1);
    end

end

% Calculate the total flux (solar + emission + solar reflection) for the requested element range:
% Load the view factor matrix for the current portion of the map (elementRange):
if exist([settings.dirPath.output, 'VFM/VFM', num2str(elementRange(1)), '.mat'], 'file')	
	load([settings.dirPath.output, 'VFM/VFM', num2str(elementRange(1)), '.mat']);
    
    % Loading the reflection and emission matrices as matfiles for faster memory
    % handling:
    tempReflectionMatrix = fRef.reflectionFlux(:,:,timeStep);
    tempEmissionMatrix = fEmis.emissionFlux(:,:,timeStep);
    
	for facetLinearIndex = 1:length(elementRange)
		tVFM = reshape(VFM(:,facetLinearIndex),mapProp.mapSize,mapProp.mapSize);
        
		% The radiation reaching the examined grid facet, originated in all other facets of the grid (i.e., both reflection and emission).
        % Visual scattering:
        reflectionReachingExaminedFacet = physProp.albedo .* (sum(sum((solarFluxMatrix + tempReflectionMatrix) .* full(tVFM))));
        % IR scattering:
        emissionReachingExaminedFacet = sum(sum(((1 - physProp.meanEmissivity) .* tempEmissionMatrix + constants.sb .* physProp.meanEmissivity .* (TsurfEqui.^4)).* full(tVFM)));
			
		% The total flux reaching the grid facet is:
		buffTotalFlux(facetLinearIndex) = ( 1 - physProp.albedo ) .* (solarFluxMatrix(elementRange(facetLinearIndex)) + reflectionReachingExaminedFacet) + physProp.meanEmissivity .* emissionReachingExaminedFacet;
        
        % Saving the flux to a temporary file, to take care of cluster
        % execution:
		bufEmission(facetLinearIndex) = emissionReachingExaminedFacet;
		bufVisReflection(facetLinearIndex) = reflectionReachingExaminedFacet;
	end
else
		% If no VFM exist for the above range, save NaNs instead:
        buffTotalFlux = NaN(size(buffTotalFlux));
end

% Saving to a temporary file to be consolidated later:
save([settings.dirPath.output, 'buffTotalFlux/buffTotalFlux_', num2str(elementRange(1)),'.mat'], 'buffTotalFlux');
save([settings.dirPath.output, 'buffTotalFlux/bufEmission_', num2str(elementRange(1)),'.mat'], 'bufEmission');
save([settings.dirPath.output, 'buffTotalFlux/bufVisReflection_', num2str(elementRange(1)),'.mat'], 'bufVisReflection');

% Open the lock file and append the char "1":
incrementLock(simDir);
