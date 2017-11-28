function[shadowMatrix, shadowDepthMatrix, solarFluxMatrix] = calcShadowMatrix(solarIncidenceAngle, solarAzimuth, simDir)
% This function calculates BOTH the shadow function and the resulting
% insolation map on the surface.

% Load configuration files and assign variables:
load([simDir, '/config/bin/settings.mat']);
load([settings.dirPath.config, 'bin/constants.mat']);
load([settings.dirPath.config, 'bin/physProp.mat']);
load([settings.dirPath.config, 'bin/mapProp.mat']);

% Quit the script if there is already a solar flux matrix for the current timestep:
if exist([settings.dirPath.output, 'Shadow/solarFluxMatrix_', num2str(solarAzimuth), '_', num2str(solarIncidenceAngle) ,'.mat'], 'file')
    writeToLog(['Solar flux map for incidence ', num2str(solarIncidenceAngle), ' and azimuth ', num2str(solarAzimuth), ' already exists.'], true);
    return;
end

% Loading the topography:
load([settings.dirPath.input,'Z.mat']);

% Initialize the shadow matrix:
shadowMatrix = zeros(size(Z));
shadowDepthMatrix = zeros(size(Z));

% Else, continue calculating the shadow matrix:

[X,Y] = meshgrid(-fix((length(Z) / 2)):fix((length(Z) / 2)), -fix((length(Z) / 2)):fix((length(Z) / 2)));
minX = min(X(:));
maxX = max(X(:));

% Slopes and aspect:
[dX, dY] = gradient(Z);
[dTHETA, dR] = cart2pol(dX, dY);
slopeField = atand(dR);
slopeAspectField = radtodeg(dTHETA) - 180;

%% More calculations:
azimuthAngleQuadrant = mod(floor([solarAzimuth]/90), 4) + 1;

%%
%
%  Beginning of the shadowing model:
%

%%  Calculating the solar shading matrix:
writeToLog(['Calculating shadow matrix for Zenith angle ', num2str(solarIncidenceAngle), ' and azimuth ', num2str(solarAzimuth),'.']);
% Iterate over all the facets in the Z grid, and check which
% one is hidden from the sun:
for facetLinearIndex = 1:numel(Z)
    x = @(y) X(facetLinearIndex) + (y - Y(facetLinearIndex))./tand(solarAzimuth);
    y = @(x) Y(facetLinearIndex) + tand(solarAzimuth) * (x - X(facetLinearIndex));
    
    if (solarAzimuth > 0 && solarAzimuth < 45)    
        idx = Y == floor(y(X));
        
    elseif (solarAzimuth > 45 && solarAzimuth <= 90)
        idx = X == floor(x(Y));
        
    elseif (solarAzimuth > 90 && solarAzimuth <= 135)
        idx = X == ceil(x(Y));
        
    elseif (solarAzimuth > 135 && solarAzimuth <= 180)
        idx = Y == floor(y(X));
        
	elseif (solarAzimuth > 180 && solarAzimuth <= 225)
        idx = Y == ceil(y(X));
        
	elseif (solarAzimuth > 225 && solarAzimuth <= 270)
        idx = X == ceil(x(Y));
        
	elseif (solarAzimuth > 270 && solarAzimuth <= 315)
        idx = X == floor(x(Y));
        
	elseif (solarAzimuth > 315 && solarAzimuth <= 360)
        idx = Y == floor(y(X));
    end
    
    if (azimuthAngleQuadrant == 1)
        idx(X < X(facetLinearIndex) | Y < Y(facetLinearIndex)) = 0;
    elseif (azimuthAngleQuadrant == 2)
        idx(X > X(facetLinearIndex) | Y < Y(facetLinearIndex)) = 0;
    elseif (azimuthAngleQuadrant == 3)
        idx(X > X(facetLinearIndex) | Y > Y(facetLinearIndex)) = 0;
    elseif (azimuthAngleQuadrant == 4)
        idx(X < X(facetLinearIndex) | Y > Y(facetLinearIndex)) = 0;
    end
    
    % zVector is the vector of Z(X,Y) on the line connecting the examined
    % facet to the sun;
    zVector = Z(idx);
    
    % Finally, calculating the slope:
    heightDiffVector = zVector - Z(facetLinearIndex);
    distanceVector = sqrt((X(idx) - X(facetLinearIndex)).^2 + (Y(idx) - Y(facetLinearIndex)).^2);
    slopeTangentToOtherPoint = heightDiffVector ./ distanceVector;
    
    % Compare the "slope tangent to other point" vector to the tangent of
    % the incidence angle:
    logicalSlopeComparison = slopeTangentToOtherPoint(1:end) >= tand(90 - solarIncidenceAngle);
    
    if any(logicalSlopeComparison(2:end))
        shadowMatrix(facetLinearIndex) = true;
        
        % Create a vector of all facets casting a shadow on the facet
        % of interest.
        whoShadowed = logicalSlopeComparison .* zVector(1 :end);
        whoShadowed(whoShadowed == 0) = NaN;
        [maxZ, maxX] = max(whoShadowed);
        
        % Set the reference plane as the mean of the 1D profile:
        referencePlane = mean(zVector(1:end));
        shadowDepthMatrix(facetLinearIndex) = (maxZ - referencePlane) - (Z(facetLinearIndex) - referencePlane) - distanceVector(maxX) .* cotd(solarIncidenceAngle);
        true;
    else
        shadowMatrix(facetLinearIndex) = false;
        shadowDepthMatrix(facetLinearIndex) = 0;
    end
end

%% Solar flux matrix:
%  Creating the matrix of angle of incidence:
cosIncidenceMatrix = cosd(solarIncidenceAngle) .* cosd(slopeField) + sind(solarIncidenceAngle).* sind(slopeField) .* cosd(solarAzimuth - slopeAspectField);

solarFluxMatrix = (constants.earthSolarConstant ./ settings.finiteSunArea )./ physProp.distanceFromSun.^2 .* cosIncidenceMatrix .* ~(shadowMatrix);

%  Removing values lower than zero (due to negative incidence angle).
solarFluxMatrix(solarFluxMatrix < 0) = 0;
return;

