function randomize_rx_positions( h_layout, max_dist, min_height, max_height, track_length, rx_ind )
%RANDOMIZE_RX_POSITIONS Generates random Rx positions and tracks
%
%   RANDOMIZE_RX_POSITIONS( max_dist , min_height , max_height , track_length) 
%   Places the users in the layout at random positions. Each user will be
%   assigned a linear track with random direction. The random height of the
%   user terminal will be in between 'min_height' and 'max_height'.
%
% The input variables are:
%	"max_dist":
%       the maximum distance from the layout center in [m]. Default is 50 m. 
%
%	"min_height":
%       the minimum user height in [m]. Default is 1.5 m.
%
%	"max_height":
%       the maximum user height in [m]. Default is 1.5 m.
%
%   "track_length":
%       the length of the linear track in [m]. Default is 1 m.
%
%
% QuaDRiGa Copyright (C) 2011-2014 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
% 
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published 
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% Parse input variables
if ~exist( 'max_dist' , 'var' ) || isempty( max_dist )
    max_dist = 50;
elseif ~( all(size(max_dist) == [1 1]) &&...
        isnumeric(max_dist) &&...
        isreal(max_dist) &&...
        max_dist > 0 )
    error('??? "max_dist" must be a real scalar  > 0')
end

if ~exist( 'min_height' , 'var' ) || isempty( min_height )
    min_height = 1.5;
elseif ~( all(size(min_height) == [1 1]) &&...
        isnumeric(min_height) &&...
        isreal(min_height) )
    error('??? "min_height" must be a real scalar')
end

if ~exist( 'max_height' , 'var' ) || isempty( max_height )
    max_height = 1.5;
elseif ~( all(size(max_height) == [1 1]) &&...
        isnumeric(max_height) &&...
        isreal(max_height) )
    error('??? "max_height" must be a real scalar  > 0')
end

if ~exist( 'track_length' , 'var' ) || isempty( track_length )
    track_length = 1;
elseif ~( all(size(track_length) == [1 1]) &&...
        isnumeric(track_length) &&...
        isreal(track_length) )
    error('??? "track_length" must be a real scalar  > 0')
end

if ~exist( 'rx_ind' , 'var' ) || isempty( rx_ind )
    rx_ind = 1:1:h_layout.no_rx;
end


% Generate random positions and tracks
rx_position_new = h_layout.rx_position;
for i_rx = 1:numel( rx_ind )
    n = rx_ind( i_rx );
    
    a = (2*rand-1)*max_dist + 1j*(2*rand-1)*max_dist;
    while abs(a)>max_dist
        a = (2*rand-1)*max_dist + 1j*(2*rand-1)*max_dist;
    end
    
    b = rand * (max_height - min_height) + min_height;
    
    rx_position_new(:,n) = [ real(a) ; imag(a) ; b ];
    
    h_layout.track(n).generate( 'linear',track_length );
    h_layout.track(n).name = ['Rx',num2str(n,'%04.0f')];
    if track_length>0
        h_layout.track(n).interpolate_positions( h_layout.simpar.samples_per_meter );
        h_layout.track(n).compute_directions;
    end
    scenarios = cell(h_layout.no_tx,1);
    for m=1:h_layout.no_tx
        scenarios(m) = h_layout.track(n).scenario(1);
    end
    h_layout.track(n).scenario = scenarios;
end

h_layout.rx_position = rx_position_new;
end