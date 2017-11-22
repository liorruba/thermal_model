function [T1] = heatDiffusionSemiImplicit_1d(x, dt, T0, specificHeatCapacity, density, thermalConductivity, emissivity, incomingFlux, geothermalFlux)
%% Auxiliary functions:
    % Update diagonal matrices:
    function [A,B] = updateDiagonalMatrices(alpha, beta, theta, phi)
        % Building the matrices centers (the first and last row will be changed later using the BC).
        A = diag(-alpha(2:end-1), -1) +...
            diag([0; (alpha(3:end-1)+alpha(2:end-2)+beta(2:end-1)); 0], 0) +...
            diag(-alpha(2:end-1), 1);

        B = diag(alpha(2:end-1), -1) +...
            diag([0; (-alpha(3:end-1)-alpha(2:end-2)+beta(2:end-1)); 0], 0) +...
            diag(alpha(2:end-1), 1);
        
        if nargin > 2
            % First and last rows:
            % Top BC:
            A(1,1) = (alpha(2) + alpha(1) + beta(1) + alpha(1) .* theta ./ phi);
            B(1,1) = -(alpha(2) + alpha(1) - beta(1) + alpha(1) .* theta ./ phi);

            % Bottom BC:
            A(end,end) = beta(end) + alpha(end);
            B(end,end) = beta(end) - alpha(end);
        end
    end       

%% Make grids:
M = length(x) - 1;

% Make dx vector:
dx = diff(x); dx = [dx; dx(end)];

% Constants:
stefanBoltzmannConstant = 5.67e-8;
%% Temperature and Flux Grid initialization:
T1 =                zeros(M+1, 1);
inputFluxVector =   zeros(M, 1);

%% Thermal parameters grid:
% Setting density vector for layers 0 and M+1:
density = [density(1); density; density(end)];

% Setting thermal conductivity for layers 0 and M+1:
k0 = [thermalConductivity(1); thermalConductivity; thermalConductivity(end)];

% Setting the mid-point thermal conductivity vector.
kHalf = (k0(1:end-1) + k0(2:end))/2;

% Temperature dependence volumetric heat capacity:
c0 = [specificHeatCapacity(1); specificHeatCapacity; specificHeatCapacity(end)];
rho_c = density .* c0;

% Setting the mid-point thermal rho_c vector.
rho_cHalf = (rho_c(1:end-1) + rho_c(2:end))/2;

%% Before the simulation starts, populate the diagonal matrices with initial values:
% Setting the auxiliary (temperature dependent) matrix elements:
alpha = kHalf./dx;
beta = rho_cHalf(2:end) .* (dx(2:end) + dx(1:end-1)) ./ dt;

[A, B] = updateDiagonalMatrices(alpha, beta);

%% Iterate:
% Building the tri-diagonal matrix:
% This matrix will be updated once every 20 time steps, for
% computational feasibility.

% Set boundary condition parameters:
theta = 2 * emissivity * stefanBoltzmannConstant .* (T0(1)).^3 - alpha(1);
phi = 2 * emissivity * stefanBoltzmannConstant .* (T0(1)).^3 + alpha(1);
psi = 3 * emissivity * stefanBoltzmannConstant .* (T0(1)).^4;

% Top BC:
A(1,1) = (alpha(2) + alpha(1) + beta(1) + alpha(1) .* theta ./ phi); 
B(1,1) = -(alpha(2) + alpha(1) - beta(1) + alpha(1) .* theta ./ phi);

% Bottom BC:
A(end,end) = beta(end) + alpha(end);
B(end,end) = beta(end) - alpha(end);

% Setting the input flux vector:
inputFluxVector(1) = 2 .* alpha(1) .* (incomingFlux + psi) ./ phi;
inputFluxVector(end) = geothermalFlux;

% Solve the matrix equation:
T1(2:end) = A\(B*T0(2:end) + inputFluxVector);
T1(1) = (incomingFlux + psi + T1(2) .* (phi - theta)) ./ 2 ./ phi;

end