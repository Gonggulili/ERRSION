function [ sf,kf ] = get_sf_profile( h_parset, evaltrack, i_mobile )
%GET_SF_PROFILE Calculates the SF and KF along the given track
%
%   ksi =  GET_SF_PROFILE( evaltrack ) returns the shadow fading and the
%   K-factor along the given track. This function is mainly used by the channel
%   builder class to scale the output channel coefficients. The profile is
%   calculated by using the data in the correlation maps and interpolating it to
%   the positions in the given track. Increasing the resolution of the maps
%   also increases the resolution of the profile.
%
% QuaDRiGa Copyright (C) 2011-2016 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
%
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.


% Parse input variable
if ~( isa(evaltrack, 'track') )
    error('??? "evaltrack" must be of class "track".')
elseif ~any( size(evaltrack) == 1  )
    error('??? "evaltrack" must be vector.')
end

drifting_precision = h_parset.simpar.drifting_precision;
if nargin ~= 3
    drifting_precision = 0;
end

n_snapshots = evaltrack.no_snapshots;

% Get the positions along the track
if drifting_precision == 3
    gdir = evaltrack.ground_direction;
    c_gdir = cos(gdir);
    s_gdir = sin(gdir);
    
    % All Rx can have the same antenna.
    if numel( h_parset.rx_array ) == 1 && i_mobile > 1
        e_rx = h_parset.rx_array(1).element_position;
    else
        e_rx = h_parset.rx_array(i_mobile).element_position;
    end
    
    e_rx_x = e_rx(1,:);
    e_rx_y = e_rx(2,:);
    n_rx = size(e_rx,2);
    
    % Apply the rotation
    q_s = zeros( 2,n_snapshots,n_rx );
    for i_rx = 1:n_rx
        q_s(1,:,i_rx) = c_gdir.*e_rx_x(i_rx) - s_gdir.*e_rx_y(i_rx);
        q_s(2,:,i_rx) = s_gdir.*e_rx_x(i_rx) + c_gdir.*e_rx_y(i_rx);
    end
    
    x = evaltrack.positions(1,:) + evaltrack.initial_position(1);
    x = q_s(1,:) + repmat( x,[1 n_rx] );
    
    y = evaltrack.positions(2,:) + evaltrack.initial_position(2);
    y = q_s(2,:) + repmat( y,[1 n_rx] );
    
else
    x = evaltrack.positions(1,:) + evaltrack.initial_position(1);
    y = evaltrack.positions(2,:) + evaltrack.initial_position(2);
    n_rx = 1;
end

par = evaltrack.par;
if isempty( par ) || isempty(par.pg) || isempty(par.kf)
    % If there are no precalculated values for PG or KF given, interpolate
    % the maps to obtain them.
    
    % If parts of the track lie outside the map, stop with an error message
    if min(x) < h_parset.map_extent(1, 1) || ...
            max(x) > h_parset.map_extent(1, 2) || ...
            min(y) < h_parset.map_extent(2, 1) || ...
            max(y) > h_parset.map_extent(2, 2)
        error('Parts of the track are outside the map boundaries.')
    end
    
    % Get the segment of the map where the positions lie in
    mi = find(h_parset.map_x_coord < min(x), 3, 'last');
    ma = find(h_parset.map_x_coord > max(x), 3, 'first');
    px = mi(1):ma(end);
    
    % Get the segment of the map where the positions lie in
    mi = find(h_parset.map_y_coord < min(y), 3, 'last');
    ma = find(h_parset.map_y_coord > max(y), 3, 'first');
    py = mi(1):ma(end);
    
    % Extract a subset of the maps (k-factor and shadow fading)
    z = ones( numel(py), numel(px), 2 );
    if numel( h_parset.kf_map ) == 1
        z(:,:,1) = z(:,:,1) .* h_parset.kf_map;
    else
        z(:,:,1) = h_parset.kf_map( py,px );
    end
    if numel( h_parset.sf_map ) == 1
        z(:,:,2) = z(:,:,2) .* h_parset.sf_map;
    else
        z(:,:,2) = h_parset.sf_map( py,px );
    end
    
    % Use a custom 2D spline interpolation to interpolate the maps along the
    % track.
    ksi = spline_2d( [x;y] , z , h_parset.map_x_coord(px) , h_parset.map_y_coord(py) );
    
    % Format the output
    kf = 10.^( 0.1*ksi(1,:) );
    kf = reshape( kf,n_snapshots,n_rx);
    
    sf = 10.^( 0.1*ksi(2,:) );
    sf = reshape( sf,n_snapshots,n_rx);
end

% Read the SF and KF from the precalculated values in "evaltrack.par"
% Previously calculated value will be overwritten (e.g. when only partial
% initial LSPs are provided in "evaltrack.par".
if ~isempty( par )
    
    if drifting_precision == 3
        o_rx = ones(1,n_rx);
        if ~isempty(par.pg)
            sf =  10.^( 0.1*par.pg( o_rx , : ).');
        end
        if ~isempty(par.kf)
            kf =  10.^( 0.1*par.kf( o_rx , : ).');
        end
        
    else
        % By default, only one vector of SF or KF values can be given here.
        % If there are more than one in 'evaltrack.par', then it is likely
        % that they belong to different Txs. However, we have no way of
        % finding out, which Tx is meant here. This is sorted out in
        % 'layout.create_parameter_sets' instead. Hence, we only use the
        % first vector (or row) and throw a warning.
        if ~isempty(par.pg)
            sf =  10.^( 0.1*par.pg(1,:).');
            
            if size( par.pg,1 ) > 1
                warning('QuaDRiGa:parameter_set:get_sf_profile',...
                    ['Multiple path gain values found. ',...
                    'There should be only one value per subtrack.']);
            end
        end
        if ~isempty(par.kf)
            kf =  10.^( 0.1*par.kf(1,:).');
        end
        
    end
end

end
