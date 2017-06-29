function gain_dBi = normalize_gain( h_array, element, gain )
%NORMALIZE_GAIN Normalizes all patterns to their gain
%
% QuaDRiGa Copyright (C) 2011-2015 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
%
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

if ~exist('element','var') || isempty( element )
    element = 1:h_array.no_elements;
end

if ~exist('gain','var') || isempty( gain )
    gain = NaN(1,numel(element));
elseif numel(gain) == 1
    gain = ones(1,numel(element))*gain;
elseif numel(gain) ~= numel(element)
    error('The number of gain values must either be 1 or match the numebr of elements.')
end

if ~(any(size(element) == 1) && isnumeric(element) ...
        && isreal(element) && all(mod(element, 1) == 0) && all(element > 0))
    error('??? "element" must be integer and > 0')
elseif any(element > h_array.no_elements)
    error('??? "element" exceed "no_elements"')
end

[~, elev_grid] = meshgrid(h_array.azimuth_grid, h_array.elevation_grid);
gain_dBi = zeros(numel(element),1);

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
    
    if isnan( gain(n) )
        gain_lin    = sum(tmp) ./ sum(P(:).*tmp);
    else
        gain_lin    = 10.^(0.1*gain(n));
    end

    gain_dBi(n) = 10*log10(gain_lin);
    
    % Normalize Patterns by their gain
    tmp = sqrt(gain_lin./P_max);
    h_array.Fa(:, :, element(n)) = Fa .* tmp;
    h_array.Fb(:, :, element(n)) = Fb .* tmp;
    if numel( h_array.Fc ) > 1
         h_array.Fc(:, :, element(n)) = Fc .* tmp;
    end
end

end
