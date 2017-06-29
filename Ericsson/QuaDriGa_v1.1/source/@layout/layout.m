classdef layout < handle & matlab.mixin.Copyable
%LAYOUT Network layout definition class
%
% DESCRIPTION
% Objects of this class define the network layout of a simulation run. Each
% network layout has one or more transmitters and one or more receivers. Each
% transmitter and each receiver need to be equipped with an antenna array which
% is defined by the array class. In general, we assume that the transmitter is at
% a fixed position and the receiver is mobile. Thus, each receivers movement is
% described by a track.
%
% EXAMPLE
%
%    a = array('dipole');                     % Generate dipole array
%
%    l = layout;                              % Create new layout
%    l.simpar.center_frequency = 2.1e9;       % Set simulation parameters
%    l.simpar.sample_density = 8;             % Set sample density
%
%    l.no_tx = 2;                             % We want two Tx
%    l.tx_position = [-50 50 ; 0 0 ; 30 30];  % Tx are at 30m height and 100m apart
%    l.tx_array = a;                          % All Tx have a dipole antenna
%
%    l.no_rx = 10;                            % 10 Receivers
%    l.randomize_rx_positions( 300,1,2 );     % Rx radius: 300m, height: 1-2m
%    l.track.set_scenario({'C2l','C2n'});     % Assign scenarios to the Rx
%    l.rx_array = a;                          % All Rx have a dipole antenna
%
%    l.set_pairing;                           % Evaluate all links
%    [p,t,rx,tx] = l.create_parameter_sets;   % Generate input for channel_builder
%
%
% QuaDRiGa Copyright (C) 2011-2013 Fraunhofer Heinrich Hertz Institute
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
        name = 'Layout';          	% Name of the layout
        simpar                   	% Handle of a 'simulation_parameters' object
    end
    
    properties(Dependent)
        no_tx                       % Number of transmitters (or base stations)
        no_rx                     	% Number of receivers (or mobile terminals)
        tx_name                    	% Identifier of each Tx, must be unique
        tx_position               	% Position of each Tx in global Cartesian coordinates using units of [m]
        tx_array                  	% Handles of 'array' objects for each Tx
        rx_name                   	% Identifier of each Tx, must be unique
        rx_position               	% Initial position of each Rx (relative to track start)
        rx_array                   	% Handles of array objects for each Rx
        track                    	% Handles of track objects for each Rx
        
        % An index-list of links for which channel are created. The first
        % row corresponds to the Tx and the second row to the Rx. 
        pairing
    end
    
    properties(Dependent,SetAccess=protected)
        % Number of links for which channel coefficients are created (read only)
        no_links                    
    end
    
    properties(Access=private)
        Pno_tx          = 1;
        Pno_rx          = 1;
        Ptx_name        = {'Tx01'};
        Ptx_position    = [0;0;25];
        Ptx_array
        Prx_array
        Ptrack
        Ppairing        = [1;1];
    end
    
    properties(Hidden)
        trans_local_wgs84 = [];
        trans_wgs84_local = [];
    end
    
    methods(Access = protected)
        function cpObj = copyElement(obj)
            cpObj = copyElement@matlab.mixin.Copyable(obj);
            % deep copy of all handle objects within the class properties
            cpObj.simpar = obj.simpar.copy;
            cpObj.tx_array = obj.tx_array.copy_objects;
            cpObj.rx_array = obj.rx_array.copy_objects;
            cpObj.track = obj.track.copy_objects;
        end
    end
    
    methods
        % Constructor
        function h_layout = layout( simpar )
            
            % At some times, MATLAB seems to use old values from the memory
            % when constructing a new layout. We prevent this by assigning
            % the default settings here.
            h_layout.name = 'Layout';
            h_layout.Pno_tx = 1;
            h_layout.Pno_rx = 1;
            h_layout.Ptx_name = {'Tx01'};
            h_layout.Ptx_position = [0;0;25];
            h_layout.Ptrack = track('linear');
            h_layout.rx_name = {'Rx1'};
            h_layout.Ptx_array = array('omni');
            h_layout.Prx_array = array('omni');

            if nargin >= 1
                h_layout.simpar = simpar;
            else
                h_layout.simpar = simulation_parameters;
            end
        end
        
        % Get functions
        function out = get.no_tx(obj)
            out = obj.Pno_tx;
        end
        function out = get.no_rx(obj)
            out = obj.Pno_rx;
        end
        function out = get.tx_name(obj)
            out = obj.Ptx_name;
        end
        function out = get.tx_position(obj)
            out = obj.Ptx_position;
        end
        function out = get.tx_array(obj)
            out = obj.Ptx_array;
        end
        function out = get.rx_name(obj)
            out = cat( 2 , {obj.track.name} );
        end
        function out = get.rx_position(obj)
            out = cat( 2, obj.track.initial_position );
        end
        function out = get.rx_array(obj)
            out = obj.Prx_array;
        end
        function out = get.track(obj)
            out = obj.Ptrack;
        end
        function out = get.pairing(obj)
            out = obj.Ppairing;
        end
        function out = get.no_links(obj)
            out = size( obj.pairing , 2 );
        end
        
        % Set functions
        function set.name(obj,value)
            if ~( ischar(value) )
                error('??? "name" must be a string.')
            end
            obj.name = value;
        end
        
        function set.simpar(obj,value)
            if ~( isa(value, 'simulation_parameters') )
                error('??? "simpar" must be objects of the class "simulation_parameters".')
            elseif ~all( size(value) == [1,1]  )
                error('??? "simpar" must be scalar.')
            end
            obj.simpar = value;
        end
        
        function set.no_tx(obj,value)
            if ~( all(size(value) == [1 1]) && isnumeric(value) ...
                    && isreal(value) && mod(value,1)==0 && value > 0 )
                error('??? "no_tx" must be integer and > 0')
            end
            
            if obj.no_tx > value
                obj.Ptx_name                    = obj.Ptx_name(1:value);
                obj.Ptx_position                = obj.Ptx_position(:,1:value);
                obj.Ptx_array                   = obj.Ptx_array(1:value);
                
                ind = obj.pairing(1,:)<=value;
                obj.pairing = obj.pairing( :,ind );
                
            elseif obj.no_tx < value
                new_name = cell( 1 , value );
                for n = 1:value
                    if n <= obj.no_tx
                        new_name{n} = obj.Ptx_name{n};
                    else
                        new_name{n} = ['Tx',num2str(n,'%02u')];
                    end
                end
                obj.Ptx_name = new_name;
                obj.Ptx_position = [ obj.Ptx_position,...
                    [ zeros( 2 , value-obj.no_tx ) ; ones( 1 , value-obj.no_tx )*25 ] ];
                for n = obj.no_tx+1 : value
                    obj.Ptx_array(n) = obj.Ptx_array( 1 );
                end
            end
            obj.Pno_tx = value;
            obj.set_pairing('all');
        end
        
        function set.no_rx(obj,value)
            if ~( all(size(value) == [1 1]) && isnumeric(value) ...
                    && isreal(value) && mod(value,1)==0 && value > 0 )
                error('??? "no_rx" must be integer and > 0')
            end
            
            if obj.no_rx > value
                obj.Prx_array                   = obj.Prx_array(1:value);
                obj.Ptrack                      = obj.Ptrack(1:value);
                
                ind = obj.pairing(2,:)<=value;
                obj.pairing = obj.pairing( :,ind );
                
            elseif obj.no_rx < value
                for n = obj.no_rx+1 : value
                    obj.Prx_array(n) = obj.Prx_array( 1 );
                    obj.Ptrack(n) = track;
                    obj.Ptrack(n).name = ['Rx',num2str(n)];
                end
            end
            obj.Pno_rx = value;
            obj.set_pairing('all');
        end
        
        function set.tx_name(obj,value)
            if ~( iscell(value) )
                error('??? "tx_name" must be a cell array.')
            elseif ~any( size(value) == 1 )
                error('??? "tx_name" must be a vector on strings.')
            end
            if size(value,1)~=1
                value = value';
            end
            if size( value , 2 ) ~= obj.no_tx
                error('??? "tx_name" must match the number of Tx.')
            end
            for n = 1:obj.no_tx
                if ~ischar( value{n} )
                    error('??? Each "tx_name" must be a string.')
                end
            end
            if numel( unique( value ) ) < numel(value)
                error('??? Each "tx_name" must be unique.')
            end
            obj.Ptx_name = value;
        end
        
        function set.tx_position(obj,value)
            if ~( isnumeric(value) && isreal(value) )
                error('??? "tx_position" must consist of real numbers')
            elseif ~all( size(value,1) == 3 )
                error('??? "tx_position" must have 3 rows')
            end
            if size(value,2) ~= obj.no_tx
                obj.no_tx = size(value,2);
            end
            obj.Ptx_position = value;
        end
        
        function set.tx_array(obj,value)
            values = size(value,2);
            if ~( isa(value, 'array') )
                error('??? "tx_array" must be objects of the class array')
            elseif ~( values == obj.no_tx || values == 1 )
                error('??? "tx_array" must match "no_tx". Try to set "no_tx" first.')
            end
            
            % Check if polarization-basis is polar-spheric
            for n = 1 : values
                if ~strcmp( value( n ).polarization_basis, 'polar-spheric' )
                    value( n ) = value( n ).copy;
                    value( n ).change_pol_basis('polar-spheric');
                end
            end
            
            if values == 1
                obj.Ptx_array(1:obj.no_tx) = value;
            else
                obj.Ptx_array = value;
            end
        end
        
        function set.rx_name(obj,value)
            if ~( iscell(value) )
                error('??? "rx_name" must be a cell array.')
            elseif ~any( size(value) == 1 )
                error('??? "rx_name" must be a vector on strings.')
            end
            if size(value,1)~=1
                value = value';
            end
            if size( value , 2 ) ~= obj.no_rx
                error('??? "rx_name" must match the number of Rx.')
            end
            for n = 1:obj.no_rx
                if ~ischar( value{n} )
                    error('??? Each "rx_name" must be a string.')
                end
            end
            if numel( unique( value ) ) < numel(value)
                error('??? Each "rx_name" must be unique.')
            end
            for n = 1:size(value,2)
                obj.track(n).name = value{n};
            end
        end
        
        function set.rx_position(obj,value)
            if ~( isnumeric(value) && isreal(value) )
                error('??? "rx_position" must consist of real numbers')
            elseif ~all( size(value,1) == 3 )
                error('??? "rx_position" must have 3 rows')
            end
            no_pos = size(value,2);
            if no_pos ~= obj.no_rx
                obj.no_rx = no_pos;
            end
            
            for n = 1:no_pos
                obj.track(n).initial_position = value(:,n);
            end
        end
        
        function set.rx_array(obj,value)
            values = size(value,2);
            if ~( isa(value, 'array') )
                error('??? "rx_array" must be objects of the class array')
            elseif ~( values == obj.no_rx || values == 1 )
                error('??? "rx_array" must match "no_rx". Try to set "no_rx" first.')
            end
            
            % Check if polarization-basis is polar-spheric
            for n = 1 : values
                if ~strcmp( value( values ).polarization_basis, 'polar-spheric' )
                    value( values ) = value( values ).copy_objects;
                    value( values ).change_pol_basis('polar-spheric');
                end
            end
            
            if values == 1
                obj.Prx_array(1:obj.no_rx) = value;
            else
                obj.Prx_array = value;
            end
        end
        
        function set.track(obj,value)
            if ~( isa(value, 'track') )
                error('??? "track" must be objects of the class track')
            end
            
            if numel(value) ~= obj.no_rx
                obj.no_rx = numel(value);
            end
            
            nm = cell( 1,numel(value) );
            for n=1:numel(value)
                nm{n} = value(n).name;
            end
            
            if  numel( unique(nm) ) < numel(value)
                error('??? Track name must be unique.')
            end
            
            obj.Ptrack = value;
        end
        
        function set.pairing(obj,value)
            value_list = reshape( value,1,[] );
            if ~( isnumeric(value) &&...
                    isreal(value) &&...
                    all( mod(value_list,1)==0 ) &&...
                    size(value,1) == 2 &&...
                    all( value_list > 0 ) )
                error('??? "pairing" must be a positive integer matrix with two rows')
            elseif any( value(1,:)>obj.no_tx )
                error('??? "pairing" refers to non-existing Tx')
            elseif any( value(2,:)>obj.no_rx )
                error('??? "pairing" refers to non-existing Rx')
            end
            
            value_new = unique( value(1,:) + 1j * value(2,:) );
            if numel( value_new ) < size(value,2)
                value = [ real( value_new ) ; imag( value_new ) ];
                warning('MATLAB:layout:multiple_pairs','removed multiple entires from "pairing".');
            end
            
            obj.Ppairing = value;
        end
    end
    
    methods(Static)
        obj = generate( varargin )
        [ h_layout , trans_local_wgs84 , trans_wgs84_local ] = ...
            kml2layout( fn , simpar , trans_wgs84_local )
    end
end
