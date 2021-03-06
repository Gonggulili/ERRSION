function pos = set_satellite_pos( h_layout, rx_latitude, sat_el, sat_az, sat_height, tx_no )
%SET_SATELLITE_POS Calculates the Tx position from a satellite orbit
%
% QuaDRiGas reference coordinate system is on the surface of the
% earth. In order to use QuaDRiga for satellite links, the satellite
% position must be set. Normally, this position is given in azimuth and
% elevation relative to the users position. This function takes a satellite
% orbital position and calculates the corresponding transmitter coordinates.    
%
% pos = SET_SATELLITE_POS( h_layout, rx_latitude, sat_el, sat_az, sat_height, tx_no ) 
%
%   Input variables:
%   	"rx_latitude":
%           The receiver latitude coordinate on the earth surface in [deg].
%           Default is 52.5. 
%
%   	"sat_el":
%           Satellite elevation in [deg]. Default is 31.6. 
%
%   	"sat_az":
%           Satellite azimuth in [deg] given in compass coordinates.
%           Default is 180 (south). 
%
%   	"sat_height":
%           Satellite height in [km] relative to earth surface. 
%           Default is 35786 (GEO orbit). 
%
%   	"tx_no":
%           The tx_no in the layout for which the position should be set.
%           Default is 1. 
%
%
% QuaDRiGa Copyright (C) 2011-2013 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
%
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.


if ~exist( 'rx_latitude' , 'var' ) || isempty( rx_latitude )
    rx_latitude = 52.5;     % Berlin :-)
elseif ~( isnumeric(rx_latitude) &&...
        isreal(rx_latitude) &&...
        all( size(rx_latitude) == 1 ) &&...
        rx_latitude >= -90 && rx_latitude <= 90 )
    error('??? "rx_latitude" has wrong format.');
end
rx_latitude = abs( rx_latitude );   % Northern and southern hemisphere are symmetric

if ~exist( 'sat_el' , 'var' ) || isempty( rx_latitude )
    sat_el = 31.6;
elseif ~( isnumeric(sat_el) &&...
        isreal(sat_el) &&...
        all( size(sat_el) == 1 ) &&...
        sat_el >= 0 && sat_el <= 90 )
    error('??? "sat_el" has wrong format.');
end

if ~exist( 'sat_az' , 'var' ) || isempty( sat_az )
    sat_az = 180;           % South
elseif ~( isnumeric(sat_az) &&...
        isreal(sat_az) &&...
        all( size(sat_az) == 1 ) &&...
        sat_az >= 0 && sat_az <= 360 )
    error('??? "sat_az" has wrong format.');
end

if ~exist( 'sat_height' , 'var' ) || isempty( sat_height )
    sat_height = 35786;     % GEO
elseif ~( isnumeric(sat_height) &&...
        isreal(sat_height) &&...
        all( size(sat_height) == 1 ) &&...
        sat_height >= 0 )
    error('??? "sat_height" has wrong format.');
end

if ~exist( 'tx_no' , 'var' ) || isempty( tx_no )
    tx_no = 1;
elseif ~( isnumeric(tx_no) &&...
        isreal(tx_no) &&...
        all( size(tx_no) == 1 ) &&...
        sum( tx_no == (1:h_layout.no_tx) ) == 1 )
    error('??? "tx_no" has wrong format.');
end

dist_x      = sat_height + rx_latitude/90 * 6378;       % [km]
dist_y      = (1-rx_latitude/90) * 6384;                % [km]

sat_dist    = sqrt(dist_x^2 + dist_y^2);                % [km]
sat_dist    = sat_dist*1e3;                             % [m]

sat_x = sat_dist * cosd(sat_el) * cosd( -sat_az+90 );
sat_y = sat_dist * cosd(sat_el) * sind( -sat_az+90 );
sat_z = sat_dist * sind(sat_el);

pos = [ sat_x ; sat_y ; sat_z ];

h_layout.tx_position(:,tx_no) = pos;

end
