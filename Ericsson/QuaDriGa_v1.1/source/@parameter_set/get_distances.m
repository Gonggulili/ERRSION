function dist = get_distances( obj )
%GET_DISTANCES Calculates the distances between Rx and Tx
%
%   GET_DISTANCES returns the distances between receiver (Rx) and transmitter
%   (Tx). Here, some of the core functions of "layout2link" from the
%   original WINNER code are implemented.
%
% QuaDRiGa Copyright (C) 2011-2012 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
% 
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published 
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

dist = sqrt( (obj.positions(1,:) - obj.tx_position(1) ).^2 +...
    ( obj.positions(2,:) - obj.tx_position(2) ).^2 );

end

