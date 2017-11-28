function [solarFlux] = loadSolarFlux(solarZenith, solarAzimuth)

for ii = 1:length(solarZenith)
    if solarZenith(ii) >= 90
        continue;
    end
    if exist(['output/Shadow/solarFluxMatrix_',num2str(solarAzimuth(ii)),'_',num2str(solarZenith(ii)),'.mat'])
        load(['output/Shadow/solarFluxMatrix_',num2str(solarAzimuth(ii)),'_',num2str(solarZenith(ii)),'.mat']);
    else
        solarFluxMatrix = zeros(101);
    end
    solarFlux(:,:,ii) = solarFluxMatrix;
%     ii
end
