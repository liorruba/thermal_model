function  [Z] = sphericalCrater(R, depthToDiameterRatio, gridSize)
%
% A copy of the peaks function that creates a spherical crater.
% R - the crater radius.
% depthToDiameterRatio - the depth to diameter ratio of the crater.
% gridSize - the size of the grid's length.
% center - a vector [a,b] where a and b are the horizontal coordinates of the crater
% center.
%
% Example: Z = sphericalCrater(80, 0.2, 201);
%


[X,Y] = meshgrid(-fix(gridSize/2):fix(gridSize/2),fix(gridSize/2):-1:-fix(gridSize/2));

height = depthToDiameterRatio * 2 * R;
sphereRadius = (R.^2 + height.^2) ./ (2*height);

if (height > sphereRadius)
    error('Height cannot be bigger than the sphere radius');
end
 
ind = find(sqrt(X.^2 + Y.^2) <= R);
Z = zeros(gridSize);
Z(ind) = sphereRadius - height - sqrt(sphereRadius.^2 - X(ind).^2 - Y(ind).^2);
