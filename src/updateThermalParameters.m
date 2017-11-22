function[] = updateThermalParameters(settings, T0)
%
% Update the thermal parameters (conductivity, rho_c).
%
load([settings.dirPath.config, 'bin/physProp.mat']);
load([settings.dirPath.config, 'bin/mapProp.mat']);

% Specific heat capacity:
physProp.specificHeatCapacity = -3.6125 + 2.7431.*T0 + 2.3616e-3 *T0.^2 -1.2340e-5.*T0.^3 + 8.9093e-9 .* T0.^4; % Hayne et. al, unpublished as of 2016;
if isrow(physProp.specificHeatCapacity )
    physProp.specificHeatCapacity = physProp.specificHeatCapacity'; 
end

% Density:
physProp.density = physProp.depthDensity - (physProp.depthDensity-physProp.surfaceDensity) .* exp(-mapProp.zAxis(2:end)./physProp.Hparameter); % kg m^-3, (Vasavada 2012).
if isrow(physProp.density)
    physProp.density = physProp.density'; 
end

% Thermal conductivity:
contactConductivity =   physProp.deepContactConductivity - (physProp.deepContactConductivity - physProp.surfaceContactConductivity) *...
                        (physProp.depthDensity - physProp.density) ./ (physProp.depthDensity - physProp.surfaceDensity);
physProp.thermalConductivity = contactConductivity .* (1 + physProp.radiativeConductivityFactor .* T0.^3); % W m^-1 K^-1, (Whipple, 1950).

save([settings.dirPath.config,'bin/physProp.mat'], 'physProp');
end
