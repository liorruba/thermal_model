function [Z,RMSslope_out] = gaussianRandomSurface(gridSize, RMSslope)
% Generate a Gaussian Random Field given the field's standard deviation
% ZStdv.
%
% Input:
% gridSize -    The size of the field Z(X,Y). Either a scalar to
%               produce a square matrix, or a 1x2 (alt. 2x1) vector
%               to produce a vector.
% RMSslope -    The requested RMS slope (angle).
%
%
% Output:
% Z        -    The generated field.
% RMSslope -    The generated field's RMS slope.
%
% Example:
% [Z,RMSslope_out] = createRmsSurface(201, 10)
%
% Written by Lior Rubanenko, Weizmann Institute of Science, May 2015


% The exponent of the power spectrum:
exponent = 3.9;
ZStdv = 1;

% Shuffle the random number generator to avoid repeating random numbers at startup:
rng('shuffle');

% Square grid
if (isscalar(gridSize))
    %   Create a gaussian white noise with a weighted standard deviation:
    noise = ZStdv*gridSize*randn(gridSize);

    %   Create a power spectrum of k^(-powerSpectrum):
    [k, l] = meshgrid(fix(-gridSize/2):ceil(gridSize/2) - 1);
    kcutoff = 3/4 * length(k)/2;
    kmat = sqrt(k.^2 + l.^2);
    kmat(kmat==0) = 1;

    Z_k = fftshift(fft2(noise)) .* (kmat.^(-exponent/2)) .* exp(-kmat/kcutoff);
    % Kill all domain-size wavelengths, as they induce high variability:
    Z_k(fix(length(Z_k)/2):fix(length(Z_k)/2)+2, fix(length(Z_k)/2):fix(length(Z_k)/2)+2) = 0;

    %   The inverse fourier transform generates the gaussian field:
    Z = ifft2(ifftshift(Z_k));
    [u,v]=gradient(Z);
    bufRMSslope = rms(sqrt(u(:).^2+v(:).^2));
    
    Z = Z * tand(RMSslope)/bufRMSslope;
    
    % Measure back the RMS slope on the new grid, to check it is correct:
    [u,v]=gradient(Z);
    RMSslope_out = atand(rms(sqrt(u(:).^2+v(:).^2)));

else
    error('Output must be a square grid.');

end
