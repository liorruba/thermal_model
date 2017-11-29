function[] = createRCM(elementRange, simDir)
% Load configuration files and assign variables:
load([simDir, '/config/bin/settings.mat']);
load([settings.dirPath.config, 'bin/constants.mat']);
load([settings.dirPath.config, 'bin/physProp.mat']);
load([settings.dirPath.config, 'bin/mapProp.mat']);

% Loading the topography:
load([settings.dirPath.input,'Z.mat']);

% The ray casting step resolution.
stepResolution = 1/mapProp.mapSize; % Pixels.

% Create the grid:
[X, Y] = meshgrid(-fix(mapProp.mapSize / 2):fix(mapProp.mapSize / 2), -fix(mapProp.mapSize / 2):fix(mapProp.mapSize / 2));

%%  This script calculates the ray-casting matrix (RCM) using a ray casting algorithm:
addition = ceil(length(X) / 2);
numberOfElementsInZ = numel(Z);

%   Preallocating:
RCM = false(numberOfElementsInZ, length(elementRange));

%   Find all the rays connecting each grid tile to the examined
%   grid tile:
iteratedElements = 1:numberOfElementsInZ;

%   v is a matrix of the x,y and z components of each grid tile.
v = [X(iteratedElements); Y(iteratedElements); Z(iteratedElements)];

internalIndex = 1;
fprintf('RCM progress: ');
for tileLinearIndex = elementRange

    if (mod(tileLinearIndex,100) == 0)
        progressBar = [num2str(tileLinearIndex/elementRange(end) * 100), ' %'];
        fprintf(progressBar);
        fprintf(repmat('\b',1,length(progressBar)-1));
        
    end

    %   R(i) is a matrix of the x,y and z components of the rays connecting Z(i)
    %   to the examined grid tile.
    %   u is a matrix containing the vector components the examined
    %   grid tile, replicated so it could later be substructed from
    %   v, to create R.
    u = repmat([...
        X(tileLinearIndex);...
        Y(tileLinearIndex);...
        Z(tileLinearIndex)]...
        , [1,length(v)]);
    R = u - v;
    
    
    %%  RCM
    %   c is the coefficient that later will be used to create the
    %   tracedRay vector.
    c = 0:stepResolution:1;
    
    %   Preallocation of variables.
%     tracedRayMatrix = false(3, numberOfElementsInZ, length(c) - 2);
%     surfaceHeightMatrix = false(length(R),1);
    
    
    % Traced ray matrix:
    tracedRayX = c' * R(1,:) + repmat(v(1,:), size(c'));
    tracedRayY = c' * R(2,:) + repmat(v(2,:), size(c'));
    tracedRayZ = c' * R(3,:) + repmat(v(3,:), size(c'));
    
    surfaceHeightZ = interp2(X, Y, Z, tracedRayX, tracedRayY);
    
    % Find neighbors and include them (on a flat surface this is amended by the View Factor):
    buffRCM = reshape(any(tracedRayZ(2:length(c)-1, :) < surfaceHeightZ(2:length(c)-1, :)),mapProp.mapSize,mapProp.mapSize);
    s=size(buffRCM);
    N=length(s);
    [c1{1:N}]=ndgrid(1:3);
    c2(1:N)={2};
    offsets=sub2ind(s,c1{:}) - sub2ind(s,c2{:}); offsets(offsets < 0) = []; offsets((tileLinearIndex + offsets) >= numberOfElementsInZ) = [];
    buffRCM(tileLinearIndex + offsets) = false;
    
    % Reshape the RCM back to the big RCM matrix:
    RCM(:, internalIndex) = reshape(buffRCM, numberOfElementsInZ, 1);
    
    internalIndex = internalIndex + 1;
end

% Save the variables (both Z and RCM) for future use in VFM calculation:
save([settings.dirPath.output,'RCM/RCM',num2str(elementRange(1)),'.mat'], 'RCM','-v7.3');

% Open the lock file and append the char "1":
incrementLock(simDir);
