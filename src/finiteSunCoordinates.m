function [finiteSunZenithAngle, finiteSunAzimuth] = finiteSunCoordinates(solarZenithAngle, solarAzimuth, simDir)
%
% This function calculates the finite Sun coordinates from the point-Sun
% coordinates (treating the point-Sun coordinates as the center of the
% finite Sun).
%

load([simDir, '/config/bin/settings.mat']);
load([settings.dirPath.config, 'bin/physProp.mat']);

AU = 1.49e8; % m

sunAngularRadius = rad2deg((physProp.sunRadius) ./ (physProp.distanceFromSun * AU)); % degrees

if settings.finiteSunArea == 1
    % Do nothing
    finiteSunZenithAngle(1,1) = solarZenithAngle;
    finiteSunAzimuth(1,1) = solarAzimuth;

elseif settings.finiteSunArea == 9
    finiteSunZenithAngle(1,1) = abs(solarZenithAngle - sunAngularRadius);
    finiteSunAzimuth(1,1) = solarAzimuth - sunAngularRadius;
    
    finiteSunZenithAngle(2,1) = solarZenithAngle;
    finiteSunAzimuth(2,1) = solarAzimuth - sunAngularRadius;
    
    finiteSunZenithAngle(3,1) = abs(solarZenithAngle + sunAngularRadius);
    finiteSunAzimuth(3,1) = solarAzimuth - sunAngularRadius;
    
    finiteSunZenithAngle(1,2) = abs(solarZenithAngle - sunAngularRadius);
    finiteSunAzimuth(1,2) = solarAzimuth;
    
    finiteSunZenithAngle(2,2) = solarZenithAngle;
    finiteSunAzimuth(2,2) = solarAzimuth;
    
    finiteSunZenithAngle(3,2) = abs(solarZenithAngle + sunAngularRadius);
    finiteSunAzimuth(3,2) = solarAzimuth;
    
    finiteSunZenithAngle(1,3) = abs(solarZenithAngle - sunAngularRadius);
    finiteSunAzimuth(1,3) = solarAzimuth + sunAngularRadius;
    
    finiteSunZenithAngle(2,3) = solarZenithAngle;
    finiteSunAzimuth(2,3) = solarAzimuth + sunAngularRadius;
    
    finiteSunZenithAngle(3,3) = abs(solarZenithAngle + sunAngularRadius);
    finiteSunAzimuth(3,3) = solarAzimuth + sunAngularRadius;
    
    finiteSunAzimuth = wrapTo360(finiteSunAzimuth);
else
    error("settings.finiteSunArea must be either 1 or 9 pixels.");
end
