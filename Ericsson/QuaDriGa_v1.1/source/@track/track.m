classdef track < handle & matlab.mixin.Copyable
%TRACK Movement model class
%
% DESCRIPTION
% One feature of the channel model is the continuous evolution of wireless
% channels when the terminal moves through the environment. A track describes the
% movement of a mobile terminal. It is composed of an ordered list of positions.
% During the simulation, one snapshot is generated for each position on the
% track.
%
% Along the track, wireless reception conditions may change, e.g. when moving
% from an area with LOS to a shaded area. This behavior is described by segments.
% A segment is a subset of positions that have similar reception conditions. Each
% segment is classified by a segment index (i.e. the center position of the
% segment) and a scenario. The scenario must be one of the supported scenarios in
% class parameter_set.
%
% EXAMPLE
%
%    t = track;
%    t.generate('circular',10*pi,0);
%    t.initial_position = [5;0;0];
%    t.segment_index = [1,40,90];
%    t.scenario = { 'C2l' , 'C2n' , 'C2l' };
%    t.interpolate_positions( 100 );
%    t.compute_directions;
%    t.visualize;
%
% The first line creates a new track. The second statement parameterizes the
% track as a circle with a circumference of 10*pi m or a diameter of 10 m. The
% last argument defines the tracks start point (in degree). Zero means, that it
% start at the positive x-axis. The next statement (t.segment_index) defines
% three segments, one that starts at the beginning, one that starts at point 40
% and one at point 90. Note, that each circular track has initially 129 points.
% Each segment is then given a scenario type. The last three lines interpolate
% the positions to 100 values per meter, adds the direction information and
% create a plot showing the track.
%
% Note: The value of initial position will be added to the values of positions.
% However, the correlated large scale parameters will only be calculated for the
% initial position. It is therefore recommended to set the starting position of
% the track to [0;0;0] and use the initial_position property to place the
% terminal relative to the transmitter.
%
%
%
% QuaDRiGa Copyright (C) 2011-2014 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
%
% Fraunhofer Heinrich Hertz Institute
% Wireless Communication and Networks
% Einsteinufer 37, 10587 Berlin, Germany
%
% This file is part of QuaDRiGa.
%
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% QuaDRiGa is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Lesser General Public License for more details.
%
% You should have received a copy of the GNU Lesser General Public License
% along with QuaDRiGa. If not, see <http://www.gnu.org/licenses/>.

properties
    name = 'Track';                 % Name of the track
end

properties(Dependent)
    % Position offset (will be added to positions)
    %   This position is given in global Cartesian coordinates (x,y, and
    %   z-component) in units of [m]. The initial position normally refers
    %   to the starting point of the track. If the track has only one
    %   segment, it is also the position for which the LSPs are calculated.
    %   The initial position is added to the values in the positions
    %   variable.
    initial_position
    
    no_snapshots                    % Number of positions on the track
    
    % Ordered list of position relative to the initial position
    %   QuaDRiGa calculates an instantaneous channel impulse response (also
    %   called snapshot) for each position on the track. 
    positions
    
    % Time (in sec) vs. distance (in m) for speed profile
    %   QuaDRiGa supports variable terminal speeds. This is realized by
    %   interpolating the channel coefficients at the output of the model.
    %   The variable track.movement_profile describes the
    %   movement along the track by associating a time-point with a
    %   distance-point on the track.
    movement_profile
    
    no_segments                     % Number of segments or states along the track
    segment_index                   % Starting point of each segment given as index in the positions vector
    
    % Scenarios for each segment along the track
    %   This variable contains the scenario names for each segment as a
    %   cell array of strings. A list of supported scenarios can be
    %   obtained by calling "parameter_set.supported_scenarios". If there
    %   is only one transmitter (i.e. one base station), the cell array has
    %   the dimension [1 x no_segments]. For multiple transmitters, the
    %   rows of the array may contain different scenarios for each
    %   transmitter. For example, in a multicell setup with three
    %   terrestrial base stations, the propagation conditions may be
    %   different to all BSs. The cell arrays than has the dimension [3 x
    %   no_segments].
    scenario
    
    % Manual parameter settings
    %   This variable contains a structure with LSPs. This can be used for
    %   assigning LSPs directly to the channel builder, e.g. when they are
    %   obtained from measurements. The structure contains the following
    %   fields:    
    %
    %      ds - The delay spread in [s] per segment
    %      kf - The Rician K-Factor in [dB] per snapshot
    %      pg - The effective path gain excluding antenna gains in [dB] per snapshot
    %      asD - The azimuth angle spread in [deg] per segment at the transmitter
    %      asA - The azimuth angle spread in [deg] per segment at the receiver
    %      esD - The elevation angle spread in [deg] per segment at the transmitter
    %      esA - The elevation angle spread in [deg] per segment at the receiver
    %      xpr - The NLOS cross-polarization in [dB] per segment
    %
    %   If there is only a subset of variables (e.g. the angle spreads are
    %   missing), then the corresponding fields can be left empty. They
    %   will be completed by the parameter sets. See also the method
    %   "track.generate_parameters" on how to fill this structure automatically. 
    par                             
    
    % Azimuth orientation of the terminal antenna for each snapshot 
    %   This variable can be calculated automatically from the positions by
    %   the function "track.compute_directions".  
    ground_direction            
    
    % Elevation orientation of the terminal antenna for each snapshot 
    %   This variable can be calculated automatically from the positions by
    %   the function "track.compute_directions".  
    height_direction
end

properties(Dependent,SetAccess=protected)
    closed                          % Indicates that the track is a closed curve
end

properties(Dependent,Hidden)
    par_nocheck
end

properties(Access=private)
    Pno_snapshots       = 1;
    Pno_segments        = 1;
    Pinitial_position   = [0;0;0];
    Ppositions          = [0;0;0];
    Pmovement_profile   = [];
    Pground_direction   = [];
    Pheight_direction   = [];
    Psegment_index      = 1;
    Pscenario           = {''};
    Ppar                = [];
    Plength             = [];
end

methods
    % Constructor
    function h_track = track( track_type , varargin )
        if nargin > 0
            h_track.generate( track_type , varargin{:} );
        else
            h_track.generate( 'linear' );
        end
    end
    
    % Get functions
    function out = get.no_snapshots(obj)
        out = obj.Pno_snapshots;
    end
    function out = get.no_segments(obj)
        out = obj.Pno_segments;
    end
    function out = get.initial_position(obj)
        out = obj.Pinitial_position;
    end
    function out = get.positions(obj)
        out = obj.Ppositions;
    end
    function out = get.movement_profile(obj)
        out = obj.Pmovement_profile;
    end
    function out = get.segment_index(obj)
        out = obj.Psegment_index;
    end
    function out = get.scenario(obj)
        out = obj.Pscenario;
    end
    function out = get.par(obj)
        out = obj.Ppar;
    end
    function out = get.ground_direction(obj)
        out = obj.Pground_direction;
    end
    function out = get.height_direction(obj)
        out = obj.Pheight_direction;
    end
    function out = get.closed(obj)
        if obj.no_snapshots>1 && ...
                all( obj.positions(:,1) == obj.positions(:,obj.no_snapshots) )
            out = true;
        else
            out = false;
        end
    end
    
    
    % Set functions
    function set.name(obj,value)
        if ~( ischar(value) )
            error('QuaDRiGa:Track:wrongInputValue','??? Name must be a string.')
        end
        obj.name = value;
    end
    
    function set.no_snapshots(obj,value)
        if ~( all(size(value) == [1 1]) && isnumeric(value) ...
                && isreal(value) && mod(value,1)==0 && value > 0 )
            error('QuaDRiGa:Track:wrongInputValue','??? "no_snapshots" must be integer and > 0')
        end
        
        if obj.Pno_snapshots > value
            obj.Ppositions = obj.Ppositions(:,1:value);
            if ~isempty( obj.Pheight_direction )
                obj.Pheight_direction = obj.Pheight_direction(:,1:value);
            end
            if ~isempty( obj.Pground_direction )
                obj.Pground_direction = obj.Pground_direction(:,1:value);
            end
            
        elseif obj.Pno_snapshots < value
            no_pos_new = value-obj.Pno_snapshots;
            pos_new = [ ones(1,no_pos_new) * obj.Ppositions( 1,obj.no_snapshots ); ...
                ones(1,no_pos_new) * obj.Ppositions( 2,obj.no_snapshots );...
                ones(1,no_pos_new) * obj.Ppositions( 3,obj.no_snapshots )];
            obj.Ppositions = [ obj.Ppositions , pos_new ];
            
            obj.Pground_direction = [];
            obj.Pheight_direction = [];
            
            obj.segment_index = obj.Psegment_index( obj.Psegment_index <= obj.no_snapshots );
        end
        
        obj.Plength = [];
        obj.Ppar = [];
        obj.Pno_snapshots = value;
    end
    
    function set.no_segments(obj,value)
        if ~( all(size(value) == [1 1]) && isnumeric(value) ...
                && isreal(value) && mod(value,1)==0 && value > 0 )
            error('QuaDRiGa:Track:wrongInputValue','??? "no_segments" must be integer and > 0')
        end
        
        if obj.Pno_segments > value
            obj.Pscenario  = obj.Pscenario(:,1:value);
            obj.Psegment_index  = obj.Psegment_index(:,1:value);
        elseif obj.Pno_segments < value
            no_new_segments = value-obj.no_segments;
            new_segments = cell( size( obj.Pscenario , 1 ),no_new_segments );
            for n = 1:numel(new_segments)
                new_segments{n} = '';
            end
            obj.Pscenario = [ obj.Pscenario , new_segments ];
            obj.Psegment_index =  [ obj.Psegment_index ,...
                ones( 1,no_new_segments )*obj.Psegment_index(end) ];
        end
        
        obj.Ppar = [];
        obj.Pno_segments = value;
    end
    
    function set.initial_position(obj,value)
        if ~( isnumeric(value) && isreal(value) )
            error('QuaDRiGa:Track:wrongInputValue','??? "initial_position" must consist of real numbers')
        elseif ~all( size(value) == [3,1] )
            error('QuaDRiGa:Track:wrongInputValue','??? "initial_position" must have three rows.')
        end
        obj.Pinitial_position = value;
    end
    
    function set.positions(obj,value)
        if ~( isnumeric(value) && isreal(value) && size(value,2) > 0 )
            error('QuaDRiGa:Track:wrongInputValue','??? "positions" must consist of real numbers')
        elseif ~all( size(value,1) == 3 )
            error('QuaDRiGa:Track:wrongInputValue','??? "positions" must have three rows.')
        end
        if size( value,2 ) ~= obj.Pno_snapshots
            obj.no_snapshots = size( value,2 );
        end
        obj.Plength = [];
        obj.Ppositions = value;
    end
    
    function set.movement_profile(obj,value)
        if isempty( value )
            obj.Pmovement_profile = [];
        else
            if ~( isnumeric(value) &&...
                    isreal(value) &&...
                    size(value,1) == 2  &&...
                    all(all( value >= 0 )) )
                error('QuaDRiGa:Track:wrongInputValue','??? "movement_profile" must be a vector of real positive numbers');
            end
            if value(2,end) > obj.get_length
                error('QuaDRiGa:Track:wrongInputValue','??? maximum distance in "movement_profile" exceeds track length');
            end
            if ~issorted( value(1,:) )
                error('QuaDRiGa:Track:wrongInputValue','??? time in "movement_profile" must be sorted');
            end
            obj.Pmovement_profile = value;
            
        end
    end
    
    function set.ground_direction(obj,value)
        if ~( isnumeric(value) && isreal(value) )
            error('QuaDRiGa:Track:wrongInputValue',...
                '??? "ground_direction" must be numeric and real');
        end
        if ~isempty(value)
            if all( size(value) == [ 1 1 ] )
                value = ones(1,obj.no_snapshots) .* value;
            elseif size(value,1) == obj.no_snapshots
                value = value.';
            elseif size(value,2) ~= obj.no_snapshots
                error('QuaDRiGa:Track:wrongInputValue',...
                    '??? no. of elements in "ground_direction" does not match no. of snapshots');
            end
        end
        obj.Pground_direction = value;
        
        % This is a dummy since height_direction is not accessible from
        % the outside
        obj.Pheight_direction = zeros(1,obj.no_snapshots);
    end
    
    function set.height_direction(obj,value)
        if ~( isnumeric(value) && isreal(value) )
            error('QuaDRiGa:Track:wrongInputValue',...
                '??? "height_direction" must be numeric and real');
        end
        if ~isempty(value)
            if all( size(value) == [ 1 1 ] )
                value = ones(1,obj.no_snapshots) .* value;
            elseif size(value,1) == obj.no_snapshots
                value = value.';
            elseif size(value,2) ~= obj.no_snapshots
                error('QuaDRiGa:Track:wrongInputValue',...
                    '??? no. of elements in "height_direction" does not match no. of snapshots');
            end
        end
        obj.Pheight_direction = value;
    end
    
    function set.segment_index(obj,value)
        if ~( isnumeric(value) &&...
                isreal(value) &&...
                all( mod(value,1)==0 ) &&...
                any( size(value) == 1 ) &&...
                all( value > 0 ) )
            error('QuaDRiGa:Track:wrongInputValue','??? "segment_index" must be a vector of integer numbers having values > 0')
        elseif min(value) ~= 1
            error('QuaDRiGa:Track:wrongInputValue','??? first segment must start at index 1')
        elseif max(value) > obj.no_snapshots
            error('QuaDRiGa:Track:wrongInputValue','??? maximum segment exceeds number of snapshots')
        end
        
        % Make sure that indices are unique and sorted
        if size(value,1) ~= 1
            value = value';
        end
        value = sort( unique( value ) );
        no_segments_old = obj.no_segments;
        
        % Set the number of segments
        if numel(value) ~= obj.no_segments
            obj.no_segments = numel(value);
        end
        
        % Preserve scenarios
        if no_segments_old == 1 && numel(value) > 1
            tmp = cell( size( obj.scenario,1 ) , numel(value) );
            for n = 1:numel(value)
                tmp(:,n) = obj.scenario(:,1);
            end
            obj.Pscenario = tmp;
        end
        
        obj.Psegment_index = value;
    end
    
    function set.scenario(obj,value)
        if ~iscell(value) && ischar(value)
            value = {value};
        end
        if ~iscell(value)
            error('QuaDRiGa:Track:wrongInputValue','??? "scenario" must be a cell array.')
        end
        if size(value, 2) == 1 && obj.no_segments > 1
            for n = 2:obj.no_segments
                value(:,n) = value(:,1);
            end
        elseif size(value, 2) ~= obj.no_segments
            error('QuaDRiGa:Track:wrongInputValue','??? The number of columns in "scenario" must match the number of segments.')
        end
        for n = 1:obj.no_segments
            for m = 1:size(value, 1)
                if ~ischar(value{m, n})
                    error('QuaDRiGa:Track:wrongInputValue',...
                        '??? Each "scenario" must be a string.')
                end
            end
        end
        obj.Pscenario = value;
    end
    
    function set.par(obj,value)
        % If par is set to an empty array, use the automatic mode.
        if ~isempty( value )
            
            % Initialize
            if isempty( obj.Ppar )
                data = struct('ds',[],'kf',[],'pg',[],...
                    'asD',[],'asA',[],'esD',[],'esA',[],'xpr',[]);
            else
                data = obj.Ppar;
            end
            
            % Check if par is a struct
            if ~isstruct( value )
                error('QuaDRiGa:Track:wrongInputValue',...
                    ['??? "par" must be a structure with either fields:',...
                    ' ds, kf, pg, asD, asA, esD, esA','xpr']);
            end
            
            supported_fields = fieldnames( data );
            names = fieldnames( value );
            
            for n = 1:numel( names )
                % Check if parameter is supported
                if ~any(strcmp( names{n}, supported_fields ))
                    error('QuaDRiGa:Track:wrongInputValue',...
                        ['??? "',names{n},'" is not a valid parameter.']);
                end
                
                val = value.( names{n} );
                if isempty( val )
                    % If the parameter is empty, delete contents and exit.
                    data.( names{n} ) = [];
                    
                else
                    % Check if parameter has the right format
                    switch names{n}
                        case {'pg','kf'}
                            % PG and KF are given in [dB] per snapshot
                            if ~( isnumeric( val ) && ...
                                    all(all(~isnan( val ))) && ...
                                    size(val,2) == obj.no_snapshots && ...
                                    isreal( val ) )
                                error('QuaDRiGa:Track:wrongInputValue',...
                                    ['??? "',names{n},'" must be real and the',...
                                    ' length must match the number of s.']);
                            end
                            
                        case 'xpr'
                            % XPR is given in [dB] per segment
                            if ~( isnumeric( val ) && ...
                                    all(all(~isnan( val ))) && ...
                                    size(val,2) == obj.no_segments && ...
                                    isreal( val ) )
                                error('QuaDRiGa:Track:wrongInputValue',...
                                    ['??? "',names{n},'" must be > 0 and the',...
                                    ' length must match the number of segments.']);
                            end
                            
                        otherwise
                            % DS in in [s] and 
                            % asD, asA, esD, esA are in [deg]
                            if ~( isnumeric( val ) && ...
                                    all(all(~isnan( val ))) && ...
                                    size(val,2) == obj.no_segments && ...
                                    all(all(val>0)) )
                                error('QuaDRiGa:Track:wrongInputValue',...
                                    ['??? "',names{n},'" must be > 0 and the',...
                                    ' length must match the number of segments.']);
                            end
                            
                    end
                    % If there was no error, save the value to data.
                    data.( names{n} ) = val;
                end
            end
        else
            data = [];
        end
        obj.Ppar = data;
    end
    
    function set.par_nocheck(obj,value)
        % Faster when we know that "par" is correct
        obj.Ppar = value;
    end
    
    % Additional methods
    function pos = positions_abs( obj )
        pos = [ obj.positions(1,:) + obj.initial_position(1) ; ...
            obj.positions(2,:) + obj.initial_position(2) ;...
            obj.positions(3,:) + obj.initial_position(3)];
    end
    
end

end
