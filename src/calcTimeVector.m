function [time] = calcTimeVector(solarAzimuth, solarZenithAngle, simDir)
    
    % Load the necessary properties (distance from Sun):
    load([simDir, '/config/bin/physProp.mat']);

    % Calculate the size of the hypotenuse of the triangle created by
    % neighbouring angles in the azmith and incidence vectors (the total
    % length of the arc the sun travelled in the sky) in degrees:
    % NOTE: the spherical pyth. theorem is invalid here since the solar distance >> incidence,zimuth.
    
    totAngle = sqrt(diff(solarAzimuth).^2 + diff(solarZenithAngle).^2);
    
    % The angular velocity of the body in degrees:
    angularVelocity = 360 / physProp.rotationPeriod;
    
    time = [0 cumsum(totAngle) ./ (360/physProp.rotationPeriod)];
end