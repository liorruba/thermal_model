function[] =  createVFM(elementRange, simDir)
% Load configuration files and assign variables:
load([simDir, '/config/bin/settings.mat']);
load([settings.dirPath.config, 'bin/constants.mat']);
load([settings.dirPath.config, 'bin/physProp.mat']);
load([settings.dirPath.config, 'bin/mapProp.mat']);

% Load the surface topography Z and the Ray Casting Matrix, RCM.
load([settings.dirPath.input,'Z.mat']);
load([settings.dirPath.output,'RCM/RCM',num2str(elementRange(1)),'.mat']);

numberOfElementsInZ = numel(Z);

% X and Y matrices.
[X, Y] = meshgrid(-fix((length(Z) / 2)):fix((length(Z) / 2)), fix((length(Z) / 2)):-1:-fix((length(Z) / 2)));

% Surface normal vector matrix:
[Nx,Ny,Nz] = surfnorm(X, Y, Z);

% Surface slopes:
[u,v] = gradient(Z);
slope = atan(sqrt(u.^2+v.^2)); clear u; clear v;

% Initiallize the VFM matrix
VFM = zeros(numberOfElementsInZ, length(elementRange));

%  View factor calculation
%  Internal indexing:
internalIndex = 1;

for tileLinearIndex = elementRange
    %   Addition is the scalar that will be added to the X,Y
    %   matrices, in order to conform them with matlab's indexing
    %   method.
    addition = ceil(length(X) / 2);
    
    %   Find all the rays connecting each grid tile to the examined
    %   grid tile:
    j = 1:numberOfElementsInZ;
    
    %   v is a matrix of the x,y and z components of each grid tile.
    v = [X(j); Y(j); Z(j)];
    
    %   R(i) is a matrix of the x,y and z components of the ray connecting Z(i)
    %   to the examined grid tile.
    %   u is a matrix containing the vector components the examined
    %   grid tile, replicated so it could later be substructed from
    %   v, to create R.
    
    u = repmat([...
        X(tileLinearIndex);...
        Y(tileLinearIndex);...
        Z(tileLinearIndex)]...
        , [1,length(v)]);
    R = v - u;
    
    %   Calculate the distance to each grid tile:
    d = sqrt(sum(R.^2,1));
    
    %   The Normals matrix to the surface Z:
    N = [reshape(Nx, 1, numberOfElementsInZ); reshape(Ny, 1, numberOfElementsInZ); reshape(Nz, 1, numberOfElementsInZ)];
    
    %   A(i) is the cos(angle) between N(i) - the normal and R(i) - the rays matrix.
    cosA = sum(N .* R,1) ./ d;
    
    %   Calculating the angles between the rays and the normal to the examined
    %   tile.
    RR = zeros(size(R)) - R;
    tileNormMatrix = repmat(N(:,sub2ind(size(Z), - Y(tileLinearIndex) + addition, X(tileLinearIndex) + addition)), [1,length(v)]);
    cosAA = sum(tileNormMatrix .* RR,1) ./ d;
    
    %   Calculating the view factor matrix:
    tVFM = (cosAA .* cosA) ./ (pi * d.^2);
    tVFM(tVFM < 0) = 0;
    %   Since the result of the view factor of a tile with "itself" is NaN, these values should be zeroed out:
    tVFM(isnan(tVFM))=0;
    
    %  Rearranging inside the FULL VFM  matrix, an N^2xlength(elementRange) matrix, where each column n contains the view factors to all other N^2-1 facets:
    VFM(:,internalIndex) = reshape(tVFM, numberOfElementsInZ, 1) ./ reshape(cos(slope), numberOfElementsInZ, 1);
    internalIndex = internalIndex + 1;
end

% Multiplying VFM By RCM:
VFM = VFM.* ~RCM;

%  Save the file using the -v7.3 flag for compression (the only way to save files >2GB).
save([settings.dirPath.output,'VFM/VFM',num2str(elementRange(1)),'.mat'], 'VFM','-v7.3');

% Open the lock file and append the char "1":
incrementLock(simDir);
