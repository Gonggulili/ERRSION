function [map, x_coords, y_coords] = power_map(h_layout, scenario, usage, sample_distance, ...
    x_min, y_max, x_max, y_min, tx_power, rx_height)
%POWER_MAP Calculates a power-map for the given layout
%
%    This function calculates receive power values in [W] on a square lattice
%    at a height of 'rx_height' above the ground for the given layout. This
%    helps to predict the performance for a given setup. 
%
% The input variables are:
%	"scenario":
%       The scenario for which the map shall be created. There are four
%       options:  
%           1. A string describing the scenario. A list of supported
%              scenarios can be obtained by calling
%              'parameter_set.supported_scenarios'.
%
%           2. A cell array of strings describing the scenario for each
%              transmitter in the layout.
%
%           3. A 'parameter_set' object. This method is useful if
%              you need to edit the parameters first. For example: call 
%              'p = parameter_set('UMal')' to load the parameters. Then edit 
%              'p.scenpar' or 'p.plpar' to adjust the settings.
%
%           4. An array of 'parameter_set' objects describing the scenario
%              for each transmitter in the layout.
%       
%    "usage"
%       A string specifying the detail level. The following options are implemented:
%           'quick'     Pattern, LOS Polarization and PL (default)
%           'sf'        Pattern, LOS Polarization, PL and SF
%           'detailed'  Full Simulation (one channel per pixel)
%           'phase'     Same as quick, but the output contains the
%                       complex-valued amplitude 
%
%    "sample_distance"  Distance between sample points in [m] (default = 10 m)
%    "x_min"   x-coordinate in [m] of the top left corner
%    "y_max"   y-coordinate in [m] of the top left corner
%    "x_max"   x-coordinate in [m] of the bottom right corner
%    "y_min"   y-coordinate in [m] of the bottom right corner
%    
%    "tx_power":
%       A vector of tx-powers in [dBm] for each transmitter in the layout.
%       This power is applied to each transmit antenna in the tx-antenna
%       array. By default (if tx_power is not given), 0 dBm are assumed.
%    "rx_height":
%       Height of the receiver points in [m] (default = 0 m)
%
%
% Output variables:
%    "map"    
%       A cell array containing the power map for each tx array in the
%       layout. The power maps are given in [W] and have the dimensions 
%           [n_y_coords, n_x_coords, n_rx_elements, n_tx_elements]
%
%    "x_coords"  
%       Vector with the x-coordinates of the map in [m]
%
%    "y_coords"
%       Vector with the y-coordinates of the map in [m]
%
%
% Example:
%    h_layout = layout;
%    h_layout.simpar.center_frequency = 2.68e9;
% 
%    h_layout.tx_position(:, 1) = [0, 100, 20];
%    h_layout.tx_array(1) = array('custom', 20, 20, 0);
% 
%    h_layout.tx_position(:, 2) = [0, 100, 45];
%    h_layout.tx_array(2) = array('custom', 20, 20, 0);
%    h_layout.tx_array(2).rotate_pattern(-90, 'z');
% 
%    h_layout.rx_array = array('dipole');
% 
%    [map, x_coords, y_coords] = h_layout.power_map( 'LOSonly', 'quick', 10, ...
%        -200, 420, 1000, -420);
% 
%    figure(1)
%    imagesc(x_coords, y_coords, map{1});
%    caxis([-5, 25])
%    set(gca, 'YDir', 'normal')
% 
%    figure(2)
%    imagesc(x_coords, y_coords, map{2});
%    caxis([-5, 25])
%    set(gca, 'YDir', 'normal')
%
%
% QuaDRiGa Copyright (C) 2011-2013 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
% 
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published 
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% Default scenario
if ~exist( 'scenario' , 'var' ) || isempty(  scenario )
    scenario = h_layout.track(1).scenario;
end

% Set the usage mode
if exist( 'usage' , 'var' ) && ~isempty(  usage )
    if strcmp( usage , 'sf' ) ||...
            strcmp( usage , 'detailed' ) ||...
            strcmp( usage , 'quick' ) ||...
            strcmp( usage , 'phase' )
        % OK
    else
        error('Usage scenario not supported.')
    end
else
    usage = 'quick';
end


if ~exist( 'sample_distance' , 'var' ) || isempty( sample_distance )
    sample_distance = 10;
end 

% Check if tx-power is given
if ~exist( 'tx_power','var' ) || isempty( tx_power )
    tx_power = zeros(1, h_layout.no_tx);
elseif numel( tx_power ) == 1 && isreal( tx_power )
    tx_power = ones( 1,h_layout.no_tx ) * tx_power;
elseif any( size(tx_power) ~= [1 h_layout.no_tx] )
    error(['??? Number of columns in "tx_power" must match the number',...
        ' of transmitters in the layout.'])
elseif isnumeric( tx_power ) && isreal( tx_power )
    % OK
else
    error('??? "tx_power" has wrong format.')
end

% Check if rx_height is given
if ~exist( 'rx_height' , 'var' ) || isempty( rx_height )
    rx_height = 0;
end

if nargin <= 4
    x_min = min( h_layout.tx_position(1,:)); 
    y_max = max( h_layout.tx_position(2,:));
    x_max = max( h_layout.tx_position(1,:));
    y_min = min( h_layout.tx_position(2,:));
    
    extend = max( [ 0.33*(x_max-x_min) , 0.33*(y_max-y_min) , 200 ] );
    
    x_min = floor( (x_min - extend)/sample_distance )*sample_distance;
    x_max = ceil(  (x_max + extend)/sample_distance )*sample_distance;
    y_max = ceil(  (y_max + extend)/sample_distance )*sample_distance;
    y_min = floor( (y_min - extend)/sample_distance )*sample_distance;
end


% Get the sample grid in x and y direction
x_coords = x_min : sample_distance : x_max;
y_coords = y_max : -sample_distance : y_min;

n_x_coords = numel(x_coords);
n_y_coords = numel(y_coords);
n_coords  = n_x_coords*n_y_coords;

drifting_precision = h_layout.simpar.drifting_precision;
h_layout.simpar.drifting_precision = 0;
n_bs = h_layout.no_tx;

% Create parameter sets for each tx-position
h_parset = parameter_set.empty(0, n_bs);
check = true;

for i_bs = 1:n_bs
    if isa( scenario ,'parameter_set' ) && numel( scenario ) == 1
        h_parset(i_bs) = parameter_set;
        h_parset(i_bs).scenpar = scenario.scenpar;
        h_parset(i_bs).plpar   = scenario.plpar;
        
    elseif isa( scenario ,'parameter_set' ) && numel( scenario ) == h_layout.no_tx
        h_parset(i_bs) = parameter_set;
        h_parset(i_bs).scenpar = scenario(i_bs).scenpar;
        h_parset(i_bs).plpar   = scenario(i_bs).plpar;
        
    elseif iscell( scenario ) && numel( scenario ) == h_layout.no_tx
        try
            h_parset(i_bs) = parameter_set( scenario{i_bs},[],check );
        catch  % If shortnames are used, the above statement will fail.
            h_parset(i_bs) = parameter_set( scenario{i_bs},[],true );
        end
        
    elseif ischar( scenario )
        try
            h_parset(i_bs) = parameter_set( scenario,[],check );
        catch  % If shortnames are used, the above statement will fail.
            h_parset(i_bs) = parameter_set( scenario,[],true );
        end
        
    else
        error('Scenario definition is invalid.');
    end
    check = false;
    
    h_parset(i_bs).name = h_layout.tx_name{i_bs};
    h_parset(i_bs).simpar = h_layout.simpar;
    h_parset(i_bs).tx_position = h_layout.tx_position(:,i_bs);
    h_parset(i_bs).no_positions = n_coords;
    h_parset(i_bs).positions = [repmat(x_coords, 1, n_y_coords); ...
        reshape(repmat(y_coords, n_x_coords, 1), 1, []); ...
        rx_height .* ones(1, n_x_coords * n_y_coords)];
    
    h_parset(i_bs).samples_per_meter = max( 2/sample_distance , 1/(3*pi) );
    h_parset(i_bs).rx_track = track('linear',0,0);
    h_parset(i_bs).rx_array = h_layout.rx_array(1);
    h_parset(i_bs).tx_array = h_layout.tx_array(i_bs);
end

if strcmp( usage , 'sf' )
	% Calculate parameter maps
    h_parset.update_parameters;
    usage = 'quick';
end

map = cell.empty(0,n_bs);

switch usage
    case 'quick'
        for i_bs = 1:n_bs
            %h_channel = channel_builder.get_los_channels(h_parset(i_bs));
            coeff = channel_builder.get_los_channels(h_parset(i_bs),'single','coeff');
            
            pow = permute(abs(coeff).^2, [3, 1, 2]);
            pow = reshape(pow, n_x_coords, n_y_coords, size( coeff, 1 ), size( coeff, 2 ));
            pow = permute(pow, [2, 1, 3, 4]);
            
            % Add tx_power
            pow = pow .* 10.^( 0.1*tx_power( i_bs ) );
            
            map{i_bs} = pow;
            
        end
        
    case 'phase'
        for i_bs = 1:n_bs
            h_channel = channel_builder.get_los_channels( h_parset(i_bs));

            cf = permute( h_channel.coeff , [4,1,2,3] );
            cf = reshape( cf , n_x_coords , n_y_coords , h_channel.no_rx , h_channel.no_tx );
            cf = permute(cf,[2,1,3,4]);
            
            % Add tx_power
            cf = cf .* sqrt( 10.^( 0.1*tx_power( i_bs ) ) );
            
            map{i_bs} = cf;
        end
        
    case 'detailed'
        % Create maps and parameters for each sample point
        h_parset.update_parameters;
        
        % Calculate channels
        h_channel = h_parset.get_channels;
        
        % Calculate the maps
        for i_bs = 1:n_bs
            ind = (i_bs-1)*n_coords;
            
            pow = zeros( n_coords, h_channel(1).no_rx , h_channel(1).no_tx );
            for n = 1:n_coords
                tmp = abs(h_channel(n+ind).coeff(:,:,:,1)).^2;
                pow(n,:,:) = sum(tmp,3);
            end
            pow = reshape( pow , n_x_coords , n_y_coords , h_channel(1).no_rx , h_channel(1).no_tx );
            pow = permute(pow,[2,1,3,4]);
            
            % Add tx_power
            pow = pow .* 10.^( 0.1*tx_power( i_bs ) );
            
            map{i_bs} = pow;
        end
end

h_layout.simpar.drifting_precision = drifting_precision;

end

