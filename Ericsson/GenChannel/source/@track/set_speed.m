function set_speed( h_track , speed , no_chk )
%SET_SPEED Sets a constant speed in m/s for the entire track
%
% QuaDRiGa Copyright (C) 2011-2012 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
%
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

if nargin < 3 || ~logical( no_chk )
    if exist( 'speed','var' )
        if ~isempty( speed )
            if ~( all(size(speed) == [1 1]) && isnumeric(speed) && isreal(speed) && speed > 0 )
                error('??? Invalid sampling interval. The value must be real and > 0.')
            end
        end
    else
        speed = 1;
    end
end

if numel(h_track) > 1
    % Do for each element in array
    for n = 1:numel(h_track)
        h_track(n).set_speed( speed );
    end
else
    if isempty( speed )
        h_track.movement_profile = [];
    else
        len = h_track.get_length;
        h_track.movement_profile = [ 0 , len/speed ; 0 len ];
    end
end

end

