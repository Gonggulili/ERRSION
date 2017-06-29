function [ gain_dBi, pow_max ] = calc_gain( h_array, element  )
%CALC_GAIN Calculates the gain of the antenna array
%
%   Output:
%       gain_dBi   Normalized Gain of the antenna
%       pow_max    Maximum power in main beam direction
%
% QuaDRiGa Copyright (C) 2011-2015 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
%
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

if nargin == 1
    element = 1:h_array.no_elements;
end

if ~(any(size(element) == 1) && isnumeric(element) ...
        && isreal(element) && all(mod(element, 1) == 0) && all(element > 0))
    error('??? "element" must be integer and > 0')
elseif any(element > h_array.no_elements)
    error('??? "element" exceed "no_elements"')
end

[~, elev_grid] = meshgrid(h_array.azimuth_grid, h_array.elevation_grid);
gain_dBi = zeros(numel(element),1);
pow_max = zeros(numel(element),1);

for n = 1 : numel(element)
    
    % Read the array elements
    Fa = h_array.Fa(:, :, element(n));
    Fb = h_array.Fb(:, :, element(n));
    if numel( h_array.Fc ) == 1
        Fc = zeros( size(Fa) );
    else
        Fc = h_array.Fc(:, :, element(n));
    end
    
    % calculate radiation power pattern
    P = abs(Fa).^2 + abs(Fb).^2 + abs(Fc).^2;
    
    % Normalize by max value
    P_max = max( P(:) );
    P = P ./ P_max;
    
    % Calculate Gain
    tmp         = cos(elev_grid(:));
    gain_lin    = sum(tmp) ./ sum(P(:).*tmp);
    
    gain_dBi(n) = 10*log10(gain_lin);
    pow_max(n) = 10*log10(P_max);
end

end

