function[convFlux] = smoothOutShadow(solarFlux, timeStep, simDir)

% Load configuration files and assign variables:
load([simDir, '/config/bin/settings.mat']);
load([settings.dirPath.config, 'bin/constants.mat']);
load([settings.dirPath.config, 'bin/physProp.mat']);
load([settings.dirPath.config, 'bin/mapProp.mat']);
load([settings.dirPath.config, 'bin/settings.mat']);

% Load data files:
load([settings.dirPath.config,'bin/solarAzimuth.mat']);
load([settings.dirPath.config,'bin/solarZenithAngle.mat']);
load([settings.dirPath.config,'bin/time.mat']);

sunAngularDiameter = (physProp.sunRadius * 2) ./ (physProp.distanceFromSun * 1.49e8); % units: km/km
tau = sunAngularDiameter / ((2*pi/physProp.rotationPeriod) * cosd(mapProp.latitude));

dt = time(timeStep) - time(timeStep-1);

if (tau <= dt && dt ~= 0)
    n = floor(tau / dt);
    writeToLog(['dt for time step ', num2str(timeStep),' is larger than tau, correcting the solar flux using a Gaussian mean.'], true);
    % Create a vector containing the temporal Gaussian function:
    % Time variable:a
    t = time(timeStep) - time((timeStep - n):(timeStep + n));
    % Create an exponential weight function:
    weightFun = exp(-t.^2/(2 * tau.^2)); 

    % Normalize the weight function:
    weightFun = weightFun./(sum(weightFun));  
    
    % Calculate the convoluted flux:
    convFlux = conv2(1, weightFun, reshape(solarFlux, (mapProp.mapSize).^2, settings.timeStepsPerDaySim),'same');
else
	convFlux = reshape(solarFlux,(mapProp.mapSize).^2, settings.timeStepsPerDaySim);
end
