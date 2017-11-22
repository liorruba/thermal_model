function[T1] = fixTemperatureInstability(z, dt, T0, specificHeatCapacity, density, thermalConductivity, emissivity, lastStepIncomingRadiation, incomingRadiation, geothermalFlux, N)
% 
% This function corrects temperature instabilities that may occur due to
% abrupt changes in the incoming flux by reducing the size of the flux
% increment.
%

% Checking whether a fix is necessary by limiting the maximum temperature
% increment during one time step to 10 K:

newDt = dt / N;
newIncomingRadiationInterval = (incomingRadiation - lastStepIncomingRadiation) ./ N;

for timeStep = 1:N
    incomingRadiation = lastStepIncomingRadiation + newIncomingRadiationInterval;
    
    T1 = heatDiffusionSemiImplicit_1d(z, newDt, T0, specificHeatCapacity, density, thermalConductivity, emissivity, incomingRadiation, geothermalFlux);    
    
    % If the temperature is stable:
    if abs(T1(1) - T0(1)) < 30 || N > 8
        % Do nothing.
    else
        T1 = fixTemperatureInstability(z, newDt, T0, specificHeatCapacity, density, thermalConductivity, emissivity, lastStepIncomingRadiation, incomingRadiation, geothermalFlux,2*N);
    end
    
    lastStepIncomingRadiation = incomingRadiation;
    T0 = T1;
    if timeStep == N;
        break;
    end
end   
end
