function set_speed( obj, speed_kmh, sampling_rate_s )
%SET_SPEED Sets the mobile speed
%
%   SET_SPEED( speed_kmh, sampling_rate )
%   This method can be used to automatically calculate the sample density
%   for a given mobile speed.
%
%   Input variables:
%       speed_kmh           speed in [km/h]
%       sampling_rate_s     channel update rate in [s]
%
% QuaDRiGa Copyright (C) 2011-2013 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
%
% Fraunhofer Heinrich Hertz Institute
% Wireless Communication and Networks
% Einsteinufer 37, 10587 Berlin, Germany


obj.samples_per_meter = 1/( speed_kmh/3.6 * sampling_rate_s);

end

