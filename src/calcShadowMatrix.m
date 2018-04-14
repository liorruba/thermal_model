function[shadowMatrix, shadowDepthMatrix, solarFluxMatrix] = calcShadowMatrix(solarIncidenceAngle, solarAzimuth, simDir)
% This function calculates BOTH the shadow function AND the resulting
% insolation map on the surface.

% Load configuration files and assign variables:
load([simDir, '/config/bin/settings.mat']);
load([settings.dirPath.config, 'bin/constants.mat']);
load([settings.dirPath.config, 'bin/physProp.mat']);
load([settings.dirPath.config, 'bin/mapProp.mat']);

% Loading the topography:
load([settings.dirPath.input,'Z.mat']);

% Initialize the shadow matrix:
shadowMatrix = zeros(size(Z));
shadowDepthMatrix = zeros(size(Z));

% Else, continue calculating the shadow matrix:

[X,Y] = meshgrid(-fix((length(Z) / 2)):fix((length(Z) / 2)), -fix((length(Z) / 2)):fix((length(Z) / 2)));
maxX = max(X(:));
maxY = max(Y(:));

% Slopes and aspect:
[dX, dY] = gradient(Z);
[THETA, R] = cart2pol(dX, dY);
slopeField = atand(R);
slopeAspectField = - 180 + radtodeg(THETA);

%%
%
%  Beginning of the ray casting model:
%

%%  Calculating the solar shading matrix:
writeToLog(['Calculating shadow matrix for Zenith angle ', num2str(solarIncidenceAngle), ' and azimuth ', num2str(solarAzimuth),'.']);
% Iterate over all the facets in the Z grid, and check which
% one is hidden from the sun:
for facetLinearIndex = 1:numel(Z)
    % Calculating the linear function y = tan(theta) * x + b
    % that is passing through all the grid facets in the
    % direction of the sun.
    fx = @(y) X(facetLinearIndex) + (y - Y(facetLinearIndex)) / tand(solarAzimuth);
    fy = @(x) tand(solarAzimuth)*(x - X(facetLinearIndex)) + Y(facetLinearIndex);

    if (solarAzimuth > 45 && solarAzimuth <= 135)
        y = linspace(Y(facetLinearIndex), maxY, 101);
        x = fx(y);
        z = qinterp2(X, Y, Z, x, y);
        
    elseif (solarAzimuth > 225 && solarAzimuth <= 315)
        y = linspace(Y(facetLinearIndex), -maxY, 101);
        x = fx(y);
        z = qinterp2(X, Y, Z, x, y);
        
    elseif (solarAzimuth > 135 && solarAzimuth <= 225)
        x = linspace(X(facetLinearIndex), -maxX, 101);
        y = fy(x);
        z = qinterp2(X, Y, Z, x, y);
        
    else 
        x = linspace(X(facetLinearIndex), maxX, 101);
        y = fy(x);
        z = qinterp2(X, Y, Z, x, y);
    end

    %   And Finally, calculating the slope:
    heightDiffVector = z - Z(facetLinearIndex);
%     clf; pcolor(X,Y,Z); hold on; plot(x(~isnan(z)), y(~isnan(z)),'r','LineWidth',2); axis equal; drawnow; 
    
    distanceVector = sqrt((y - Y(facetLinearIndex)).^2 + (x - X(facetLinearIndex)).^2);
    slopeTangentToOtherPoint = heightDiffVector ./ distanceVector;
    logicalSlopeComparison = slopeTangentToOtherPoint(1:end) >= tand(90 - solarIncidenceAngle);
    
    % Start comparing the slope tangent to other point vector from the 3rd place, since
    % the first index is the facet and the second index is the closest facet
    % in order to prevent discretization errors.
    if any(logicalSlopeComparison)
        shadowMatrix(facetLinearIndex) = true;
        
        % Create a vector of all facets casting a shadow on the facet
        % of interest.
        whoShadowed = logicalSlopeComparison .* z;
        whoShadowed(whoShadowed == 0) = NaN;
        [maxZ, maxD] = max(whoShadowed);
        
        % Set the reference plane as the mean of the 1D profile:
        z0 = mean(z(~isnan(z)));
        shadowDepthMatrix(facetLinearIndex) = (maxZ - z0) - (Z(facetLinearIndex) - z0) - distanceVector(maxD) .* cotd(solarIncidenceAngle);
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
end

