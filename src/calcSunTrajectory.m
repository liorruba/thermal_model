function [azimuth,solarZenithAngle,timevector] = calcSunTrajectory(solarDeclination, latitude, numberOfTimestepsPerDay, numberOfDays, rotationPeriod)
%
% This function returns the sun horizontal coordinates given the Solar
% declination and the latitude.
%
% Usage:
% [azimuth,solarZenithAngle,timevector] = plotSunTrajectory(solarDeclination, latitude, numberOfTimestepsPerDay, numberOfDays)
%


if numberOfDays >= 1
    dt = 1/(numberOfTimestepsPerDay);
    timeUnitVector = 0:dt:(1-dt);
    hourAngle = [];
    timevector = [];
    
    for ii=1:numberOfDays
        buff = circshift(timeUnitVector, [0 fix(numberOfTimestepsPerDay.*0.25)]);
        hourAngle = [hourAngle 360*buff];
        
        clear tbuff;
          tbuff = timeUnitVector + ii-1;
          timevector = [timevector rotationPeriod*tbuff];
    end
elseif numberOfDays < 1
    error('Number of days must be greater or equal to 1');
end

solarZenithAngle = 90 - asind(sind(latitude) .* sind(solarDeclination) + cosd(hourAngle) .* cosd(latitude) .* cosd(solarDeclination));
azimuth = 180 - atan2d(cosd(hourAngle) * sind(latitude) - tand(solarDeclination) * cosd(latitude),sind(hourAngle));
