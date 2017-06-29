function angles = get_angles( h_parset )
%GET_ANGLES Calculates the departure- and arrival angles of the LOS path between Tx and Rx
%
% GET_ANGLES returns the angles between Receiver (Rx) and Transmitter
% (Tx). Here, some of the core functions of "layout2link" from the
% original WINNER code are implemented.
%
% Output:
%	"angles":
%       The number of columns corresponds to the number of rx-positions.
%       The rows contain the four angles: 
%           Azimuth of Departure at the Tx (AoD, row 1)
%           Azimuth of Arrival at the Rx (AoA, row 2)
%           Elevation of Departure at the Tx (EoD, row 3)
%           Elevation of Arrival at the Rx (EoA, row 4)
%       
%
% QuaDRiGa Copyright (C) 2011-2012 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
% 
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published 
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% Initialize memory
angles = zeros( 4,h_parset.no_positions );

% ThetaBs 
angles(1,:) = atan2( h_parset.positions(2,:) - h_parset.tx_position(2) , ...
    h_parset.positions(1,:) - h_parset.tx_position(1) );

% ThetaMs 
angles(2,:) = pi + angles(1,:);

% EaBs
angles(3,:) = atan( ( h_parset.positions(3,:) - h_parset.tx_position(3) ) ./ ...
    sqrt( (h_parset.tx_position(1) - h_parset.positions(1,:)).^2 +...
    (h_parset.tx_position(2,:) - h_parset.positions(2,:)).^2 ) );

% When Rx and Tx are at the same position, the angle is NaN
% We set it to 0 instead.
angles(3, isnan( angles(3,:) ) ) = 0;

% EaMs
angles(4,:) = -angles(3,:);

% Mapping to (-180,180) degrees
angles = angles * 180/pi;
angles = mod( angles , 360);
angles = angles - 360*floor(angles/180);

end

