function compute_directions( obj )
%COMPUTE_DIRECTIONS Calculates ground and height orientations from positions
%
%   COMPUTE_DIRECTIONS calculates the orientations of the transmitter based
%   on the positions. If we assume that the receive antenna array is fixed
%   on a car and the car moves along the track, then the antenna turns with
%   the car when the car is changing direction. This needs to be accounted
%   for when generating the channel coefficients. This function calculates
%   the orientation based on the positions and stored the output in the
%   ground_direction and height_direction field of the track object.
%
% QuaDRiGa Copyright (C) 2011-2014 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
% 
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published 
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

if numel(obj) > 1
    % Do for each element in array
    for n=1:numel(obj)
        obj(n).compute_directions;
    end
    
else
    if obj.no_snapshots<2
        obj.Pground_direction = 0;
        obj.Pheight_direction = 0;
    end
    
    P = obj.Ppositions(:,2:end) - obj.Ppositions(:,1:end-1);
    [a, e] = cart2sph(P(1, :), P(2, :), P(3, :));
    
    n = obj.no_snapshots;
    if obj.closed
        a(n) = a(1);
        e(n) = e(1);
    else
        a(n) = a(n-1);
        e(n) = e(n-1);
    end
    
    obj.ground_direction = a;
    
    % Height-directions is not implemented in the channel_builder
    % Thus, it is always zero
    obj.height_direction = e;
    
end
end