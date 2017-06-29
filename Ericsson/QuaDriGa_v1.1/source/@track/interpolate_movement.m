function dist = interpolate_movement( h_track , si, method )
%INTERPOLATE_MOVEMENT Interpolates the movement profile to a distance vector
%
%   This function interpolates the movement profile. The distance vector at
%   the output can then be used to interpolate the channel coefficients to
%   emulate varying speeds. See also the tutorial "Applying Varying Speeds
%   (Channel Interpolation)" for more details.
%
%   Input and output variables:
%   	"si":
%           the sampling interval in [seconds]
%
%   	"method":
%           selects the interpolation method. The default is 'pchip'.
%           Supported methods are:
%               nearest - nearest neighbor interpolation
%               linear  - linear interpolation
%               spline  - piecewise cubic spline interpolation
%               pchip   - shape-preserving piecewise cubic hermite interpolation
%               cubic   - same as 'pchip'
%
%       "dist":
%           Distance of each interpolated position from the start of the
%           track in [m]
%
%
% QuaDRiGa Copyright (C) 2011-2014 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
%
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

if nargin == 3
    supported_types = {'nearest','linear','spline','pchip','cubic'};
    if ~( ischar(method) && any( strcmpi(method,supported_types)) )
        str = 'Array type not found; supported types are: ';
        no = numel(supported_types);
        for n = 1:no
            str = [str,supported_types{n}];
            if n<no
                str = [str,', '];
            end
        end
        error(str);
    end
else
    method = 'pchip';
end

if ~( all(size(si) == [1 1]) && isnumeric(si) && isreal(si) && si >= 0 )
    error('??? Invalid sampling interval. The value must be real and > 0.')
end

if isempty(h_track.movement_profile)
    h_track.set_speed;
end

mp = h_track.movement_profile;

max_time = mp(1,end);
max_dist = mp(2,end);

t = 0 : si : max_time ;

dist = interp1(mp(1,:),mp(2,:),t,method);

dist( dist<0 ) = 0;
dist( dist>max_dist ) = max_dist;


end

