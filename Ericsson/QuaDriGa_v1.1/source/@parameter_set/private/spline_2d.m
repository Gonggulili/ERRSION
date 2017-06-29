function ksi = spline_2d( positions , parameter_maps , map_x_coord , map_y_coord )
%SPLINE_2D 2D spline interpolation of the parameter maps
%
% SPLINE_2D implements a custom 2D method for interpolating the parameter
% maps. This is needed since the maps are often in a lower resolution to
% save computing time.
%
%
% QuaDRiGa Copyright (C) 2011-2013 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
%
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.


% Kernel for the spline interpolation function
K = 0.5 .* ...
    [0,  2,  0,  0; ...
    -1,  0,  1,  0; ...
     2, -5,  4, -1; ...
    -1,  3, -3,  1].';

n_pos    = size(positions,2);
n_maps   = size( parameter_maps , 3 );
n_maps4  = 4*n_maps;

o_pos    = ones( 1,n_pos );
o_maps   = ones( 1,n_maps );
o_maps4  = ones( 1,n_maps4 );

samples_per_meter = (numel(map_x_coord) - 1)/...
    (map_x_coord(end) - map_x_coord(1));

i_x = x2i(positions(1,:), map_x_coord);
i_y = x2i(positions(2,:), map_y_coord);

% Check if the x and y indices are really on the map
if any(i_x < 2) || any(i_x > numel(map_x_coord)-2)
    error('x index out of bound!')
end

if any(i_y < 2) || any(i_y > numel(map_y_coord)-2)
    error('y index out of bound!')
end

% Get the neighboring positions
px = [ i_x-1 ; i_x ; i_x+1 ; i_x+2 ];
py = [ i_y-1 ; i_y ; i_y+1 ; i_y+2 ];

% Get the offset between the actual user position and the next point
% on the map.
u = ( positions(1,:) - map_x_coord(i_x) )*samples_per_meter;
v = ( positions(2,:) - map_y_coord(i_y) )*samples_per_meter;

% Get the coefficients for the spline interpolation
uK = K * [ o_pos ; u ; u.^2 ; u.^3];
vK = K * [ o_pos ; v ; v.^2 ; v.^3];

% Interpolate values for each user position
ksi = zeros( n_maps , n_pos );
for i_pos = 1 : n_pos
    z = parameter_maps( py(:,i_pos), px(:,i_pos), : );
    z = reshape( z, 4, n_maps4 );
    
    z = sum( z.*vK( : , o_maps4*i_pos ) );
    z = reshape( z, 4, n_maps );
    
    ksi(:,i_pos) = sum( z.*uK( :, o_maps*i_pos) );
end

end

% SUB FUNCTIONS
function i_x = x2i(x, x_range)
if (any(x > x_range(end))) || any((x < x_range(1)))
    error('x value has to be in x_range!');
end
i_x = floor((x - x_range(1))/(x_range(end) - x_range(1))*(numel(x_range) - 1)) + 1;
end
