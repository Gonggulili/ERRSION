classdef parameter_set < handle & matlab.mixin.Copyable
%PARAMETER_SET Generation of correlated initial parameters
%
% DESCRIPTION
% This class implements all functions that are necessary to generate and manage
% correlated large scale parameters (LSPs). It also provides interfaces for the
% channel builder. LSPs are the shadow fading, the Rician K-Factor, the RMS delay
% spread and the four angles (elevation and azimuth at the transmitter and
% receiver). This class implements some core functions of the channel model and
% the user does normally not need to interact with it. However, if parameters such
% as the parameters from the Winner tables need to be changed, here is the place
% to do so.
%
% REFERENCE
% The main functionality was taken from the Winner channel model. "Kyösti, P.;
% Meinilä, J.; Hentilä, L. & others; {IST-4-027756 WINNER II D1.1.2 v.1.1}:
% WINNER II Channel Models; 2007".
%
% EXAMPLE
% This example manually creates the correlation parameters.
%
%    l = layout;                             % Create new layout
%    l.no_rx = 10;                           % 10 Receivers ...
%    l.randomize_rx_positions(500);          % ... with random positions in 500 radius
%
%    p = parameter_set;                      % Create new parameter set
%    p.name = '';
%    p.scenario = 'C2l';                     % Set scenario
%    p.scenpar.SF_sigma = 2;                 % Manually set parameters
%    p.scenpar.PerClusterAS_D = 1;
%    p.positions = l.rx_position;            % Set Rx positions
%    p.tx_position = l.tx_position;          % Set Tx positions
%    p.update_parameters;                    % Calculate parameters
%
%    cb = channel_builder(p);
%    c = cb.get_channels;                    % Calculate output channels
%
%
%
% QuaDRiGa Copyright (C) 2011-2016 Fraunhofer Heinrich Hertz Institute
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
        name = 'parameter_set';     % Name of the parameter set object
        simpar                      % Object of class simulation_parameters
        tx_array                    % Handles of Antenna array objects for each Tx
        rx_array                    % Handles of Antenna array objects for each Rx
        rx_track                    % Handles of Track objects for each Rx
    end
    
    % Dependent properties
    properties(Dependent)
        scenario                    % The name of the scenario (text string)
        scenpar                     % The parameter table
        plpar                       % Parameters for the path loss (scenario-dependent)
        
        % Number of receiver positions associated to this 'parameter_set' object
        %   Note that each segment in longer tracks is considered a new Rx
        %   position. 
        no_positions                
        
        % The list of initial positions for which LSPs are generated
        %   This variable is obtained from 'track.initial_position' and 'layout.rx_position' 
        positions
        
        % The transmitter position obtained from the corresponding 'layout.tx_position'
        tx_position
        
        ds                          % The RMS delay spread in [s] for each receiver position
        kf                          % The Rician K-Factor [linear scale] for each receiver position
        sf                          % The shadow fading [linear scale] for each receiver position
        asD                         % The azimuth spread of departure in [deg] for each receiver position
        asA                         % The azimuth spread of arrival in [deg] for each receiver position
        esD                         % The elevation spread of departure in [deg] for each receiver position
        esA                         % The elevation spread of arrival in [deg] for each receiver position
        xpr                         % The cross polarization ratio [linear scale] for each receiver position
       
        % The seven large-scale parameter maps in logarithmic scale
        %   Rows correspond to the y-coordinate and columns to the
        %   x-coordinate. When there is no STD of a parameter defined, then
        %   the maps will be stored as a scalar coefficient defined in the
        %   scenario parameters.
        ds_map                      % The RMS delay spread map in [log10(s)]
        kf_map                      % The Rician K-Factor [dB]
        sf_map                      % The shadow fading [dB]
        asD_map                     % The azimuth spread of departure in [log10(deg)]
        asA_map                     % The azimuth spread of arrival in [log10(deg)]
        esD_map                     % The elevation spread of departure in [log10(deg)]
        esA_map                     % The elevation spread of arrival in [log10(deg)]
        
        map_extension               % Distance in [m] that is added to each direction when generating maps
        map_extent                  % Extent of the maps in x- and y-direction [xmin, xmax; ymin, ymax] in [m] 
        map_size                    % Number of map pixels in x and y-direction [n_x_samples; n_y_samples]
        samples_per_meter           % Resolution of the LSP maps in [samples/m]
        map_valid                   % Indicates if maps contain valid data
        data_valid                  % Indicates if the data is valid
    end
    
    % The combined parameter maps
    properties(SetAccess=private,Hidden)
        parameter_maps              % The combined parameter maps
    end
    
    % Read only properties
    properties(SetAccess=private)
        LSP_xcorr_matrix            % The Cross-correlation matrix for the LSPs
        LSP_matrix_isOK             % Determines if the XCorr-matrix is positive definite
    end
    
    properties(Dependent,SetAccess=private)
        map_x_coord              	% The x-coordinates in [m] for each pixel of the maps
        map_y_coord               	% The y-coordinates in [m] for each pixel of the maps
    end
    
    % Data storage
    properties(Access=private)
        Pscenario               = '';
        Pscenpar                = [];
        Pplpar                  = [];
        Pno_positions           = 1;
        Ppositions              = [0;0;0];
        Ptx_position            = [0;0;0];
        Pds                     = 0;
        Pkf                     = 1;
        Psf                     = 1;
        PasD                    = 0;
        PasA                    = 0;
        PesD                    = 0;
        PesA                    = 0;
        Pxpr                    = 0;
        Pmap_ds                 = 0;
        Pmap_kf                 = 1;
        Pmap_sf                 = 1;
        Pmap_asD                = 0;
        Pmap_asA                = 0;
        Pmap_esD                = 0;
        Pmap_esA                = 0;
        Pmap_extension          = 20;
        Pmap_valid              = false;
        Pdata_valid             = false;
        Pmap_extent             = [-30, 30; -30, 30];
        Psamples_per_meter      = 0;
    end
    
    methods(Access = protected)
        function cpObj = copyElement(obj)
            cpObj = copyElement@matlab.mixin.Copyable(obj);
            % deep copy of selected handle objects within the class properties
            cpObj.simpar = obj.simpar.copy;
            cpObj.tx_array = obj.tx_array.copy_objects;
            cpObj.rx_array = obj.rx_array.copy_objects;
            cpObj.rx_track = obj.rx_track.copy_objects;
        end
    end
    
    methods
        
        % Constructor
        function h_parset = parameter_set( scenario, positions, check_parfiles )
            
            h_parset.simpar = simulation_parameters;
            h_parset.tx_array = array('omni');
            h_parset.rx_array = array('omni');
            h_parset.rx_track = track('linear',1,0);
            
            % Initialize Handles
            if exist( 'scenario' , 'var' )
                
                % Parse Input variables
                if exist( 'check_parfiles' , 'var' )
                    if ~( all(size(check_parfiles) == [1 1]) ...
                            && (isnumeric(check_parfiles) || islogical(check_parfiles)) ...
                            && any( check_parfiles == [0 1] ) )
                        error('??? "check_parfiles" must be 0 or 1')
                    end
                else
                    check_parfiles = true;
                end
                
                if exist( 'positions' , 'var' ) && ~isempty( positions )
                    positions_given = true;
                else
                    positions_given = false;
                end
                
                if check_parfiles    % Checking enabled
                    if positions_given
                        h_parset.positions = positions;
                        h_parset.generate_correlation_maps;
                        h_parset.update_parameters;
                    else
                        h_parset.set_scenario_table(scenario);
                    end
                else % Checking disabled
                    h_parset.set_scenario_table( scenario , false );
                    if positions_given
                        h_parset.Ppositions = positions;
                        h_parset.generate_correlation_maps;
                        h_parset.update_parameters;
                    end
                end
            end
        end
        
        % Get functions
        function out = get.scenario(obj)
            out = obj.Pscenario;
        end
        function out = get.scenpar(obj)
            out = obj.Pscenpar;
        end
        function out = get.plpar(obj)
            out = obj.Pplpar;
        end
        function out = get.no_positions(obj)
            out = obj.Pno_positions;
        end
        function out = get.positions(obj)
            out = obj.Ppositions;
        end
        function out = get.tx_position(obj)
            out = obj.Ptx_position;
        end
        
        % Get-Functions for the individual parameters
        function out = get.ds(obj)
            out = obj.Pds;
        end
        function out = get.kf(obj)
            out = obj.Pkf;
        end
        function out = get.sf(obj)
            out = obj.Psf;
        end
        function out = get.asD(obj)
            out = obj.PasD;
        end
        function out = get.asA(obj)
            out = obj.PasA;
        end
        function out = get.esD(obj)
            out = obj.PesD;
        end
        function out = get.esA(obj)
            out = obj.PesA;
        end
        function out = get.xpr(obj)
            out = obj.Pxpr;
        end

        % Get- Functions for the maps
        function out = get.ds_map(obj)
            out = obj.Pmap_ds;
        end
        function out = get.kf_map(obj)
            out = obj.Pmap_kf;
        end
        function out = get.sf_map(obj)
            out = obj.Pmap_sf;
        end
        function out = get.asD_map(obj)
            out = obj.Pmap_asD;
        end
        function out = get.asA_map(obj)
            out = obj.Pmap_asA;
        end
        function out = get.esD_map(obj)
            out = obj.Pmap_esD;
        end
        function out = get.esA_map(obj)
            out = obj.Pmap_esA;
        end
        
        % Get function for the combined maps
        function out = get.parameter_maps(obj)
            out = zeros( obj.map_size(2), obj.map_size(1), 7 );
            out(:,:,1) = obj.Pmap_ds;
            out(:,:,2) = obj.Pmap_kf;
            out(:,:,3) = obj.Pmap_sf;
            out(:,:,4) = obj.Pmap_asD;
            out(:,:,5) = obj.Pmap_asA;
            out(:,:,6) = obj.Pmap_esD;
            out(:,:,7) = obj.Pmap_esA;
        end
        
        function out = get.map_extension(obj)
            out = obj.Pmap_extension;
        end
        function out = get.map_valid(obj)
            out = obj.Pmap_valid;
        end
        function out = get.data_valid(obj)
            out = obj.Pdata_valid;
        end
        function value = get.samples_per_meter(h_parameter_set)
            value = h_parameter_set.Psamples_per_meter;
        end
        function value = get.map_x_coord(h_parameter_set)
            value = h_parameter_set.Pmap_extent(1, 1):1/h_parameter_set.samples_per_meter:h_parameter_set.Pmap_extent(1, 2);
        end
        function value = get.map_y_coord(h_parameter_set)
            value = h_parameter_set.Pmap_extent(2, 1):1/h_parameter_set.samples_per_meter:h_parameter_set.Pmap_extent(2, 2);
        end
        function value = get.map_extent(h_parameter_set)
            value = h_parameter_set.Pmap_extent;
        end
        function value = get.map_size(h_parameter_set)
            value = [numel(h_parameter_set.map_x_coord); numel(h_parameter_set.map_y_coord)];
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
        
        function set.rx_track(obj,value)
            if ~( isa(value, 'track') )
                error('??? "track" must be of class "track".')
            end
            obj.rx_track = value;
        end
        
        function set.rx_array(obj,value)
            if ~( isa(value, 'array') )
                error('??? "rx_array" must be of class "array".')
            end
            obj.rx_array = value;
        end
        
        function set.tx_array(obj,value)
            if ~( isa(value, 'array') )
                error('??? "tx_array" must be of class "array".')
            end
            obj.tx_array = value;
        end
        
        % Dependent set functions
        function set.scenario(obj,value)
            if ~( ischar(value) )
                error('??? "scenario" must be a string.')
            end
            obj.set_scenario_table(value);
        end
        
        function set.scenpar(obj,value)
            if ~( isstruct(value) )
                error('??? "scenpar" must be a structure.')
            end
            [ obj.LSP_xcorr_matrix , obj.LSP_matrix_isOK ] =...
                check_scenario_parameter_table( value );
            
            obj.map_valid = obj.map_valid && check_if_map_is_valid( obj.scenpar , value );
            obj.Pscenpar = value;
            obj.Pscenario = 'Custom';
        end
        
        function set.plpar(obj,value)
            if isstruct(value)
                obj.Pplpar = value;
            elseif isempty(value)
                obj.Pplpar = [];
            else
                error('??? "plpar" must be a structure.')
            end
        end
        
        function set.no_positions(obj,value)
            if ~( all(size(value) == [1 1]) && isnumeric(value) ...
                    && isreal(value) && mod(value,1)==0 && value >= 0 )
                error('??? "no_positions" must be integer and >= 0')
            end
            if obj.no_positions > value
                obj.Ppositions  = obj.Ppositions(:,1:value);
                obj.Pds         = obj.Pds(1:value);
                obj.Pkf         = obj.Pkf(1:value);
                obj.Psf         = obj.Psf(1:value);
                obj.PasD        = obj.PasD(1:value);
                obj.PasA        = obj.PasA(1:value);
                obj.PesD        = obj.PesD(1:value);
                obj.PesA        = obj.PesA(1:value);
                obj.Pxpr        = obj.Pxpr(1:value);
            elseif obj.no_positions < value
                val_new = value-obj.no_positions;
                obj.Ppositions  = [ obj.Ppositions  zeros( 3 , val_new ) ];
                obj.Pds         = [ obj.Pds  zeros( 1 , val_new ) ];
                obj.Pkf         = [ obj.Pkf  ones( 1 , val_new ) ];
                obj.Psf         = [ obj.Psf  ones( 1 , val_new ) ];
                obj.PasD        = [ obj.PasD  zeros( 1 , val_new ) ];
                obj.PasA        = [ obj.PasA  zeros( 1 , val_new ) ];
                obj.PesD        = [ obj.PesD  zeros( 1 , val_new ) ];
                obj.PesA        = [ obj.PesA  zeros( 1 , val_new ) ];
                obj.Pxpr        = [ obj.Pxpr  zeros( 1 , val_new ) ];
            end
            obj.Pno_positions = value;
            obj.Pdata_valid   = false;
        end
        
        function set.positions(obj,value)
            if ~( isnumeric(value) && isreal(value) )
                error('??? "positions" must consist of real numbers')
            elseif size(value,1) ~= 3
                error('??? "positions" must be have three rows')
            end
            if size(value,2) ~= obj.no_positions
                obj.no_positions = size(value,2);
            end
            obj.Ppositions = value;
        end
        
        function set.tx_position(obj,value)
            if ~( isnumeric(value) && isreal(value) )
                error('??? "tx_position" must consist of real numbers')
            elseif ~all( size(value) == [3,1] )
                error('??? "tx_position" must be have three rows and one column')
            end
            obj.Ptx_position = value;
        end
        
        function set.ds(obj,value)
            if ~( any(size(value) == 1) && isnumeric(value) ...
                    && isreal(value) && all( value >= 0 ) )
                error('??? "ds" must be real, scalar and >= 0')
            end
            if size(value,2) ~= obj.no_positions
                obj.no_positions = size(value,2);
            end
            obj.Pds = value;
        end
        
        function set.kf(obj,value)
            if ~( any(size(value) == 1) && isnumeric(value) && isreal(value) ...
                    && isreal(value) )
                error('??? "kf" must be real and scalar')
            end
            if size(value,2) ~= obj.no_positions
                obj.no_positions = size(value,2);
            end
            obj.Pkf = value;
        end
        
        function set.sf(obj,value)
            if ~( any(size(value) == 1) && isnumeric(value) ...
                    && isreal(value) )
                error('??? "sf" must be real and scalar')
            end
            if size(value,2) ~= obj.no_positions
                obj.no_positions = size(value,2);
            end
            obj.Psf = value;
        end
        
        function set.asD(obj,value)
            if ~( any(size(value) == 1) && isnumeric(value) ...
                    && isreal(value) )
                error('??? "asD" must be real and scalar')
            end
            if size(value,2) ~= obj.no_positions
                obj.no_positions = size(value,2);
            end
            obj.PasD = value;
        end
        
        function set.asA(obj,value)
            if ~( any(size(value) == 1) && isnumeric(value) ...
                    && isreal(value) )
                error('??? "asA" must be real and scalar')
            end
            if size(value,2) ~= obj.no_positions
                obj.no_positions = size(value,2);
            end
            obj.PasA = value;
        end
        
        function set.esD(obj,value)
            if ~( any(size(value) == 1) && isnumeric(value) ...
                    && isreal(value) )
                error('??? "esD" must be real and scalar')
            end
            if size(value,2) ~= obj.no_positions
                obj.no_positions = size(value,2);
            end
            obj.PesD = value;
        end
        
        function set.esA(obj,value)
            if ~( any(size(value) == 1) && isnumeric(value) ...
                    && isreal(value) )
                error('??? "esA" must be real and scalar')
            end
            if size(value,2) ~= obj.no_positions
                obj.no_positions = size(value,2);
            end
            obj.PesA = value;
        end
        
        function set.xpr(obj,value)
            if ~( any(size(value) == 1) && isnumeric(value) ...
                    && isreal(value) )
                error('??? "xpr" must be real and scalar')
            end
            if size(value,2) ~= obj.no_positions
                obj.no_positions = size(value,2);
            end
            obj.Pxpr = value;
        end
        
        function set.map_extension(obj,value)
            if ~( all(size(value) == [1 1]) && isnumeric(value) ...
                    && isreal(value) && value > 0 )
                error('??? "map_extension" must be real, scalar and > 0')
            end
            obj.Pmap_extension = value;
        end
        
        function set.ds_map(obj,value)
            if isnumeric(value) && isreal(value) && ...
                    ( all( size(value) == [obj.map_size(2), obj.map_size(1)] ) || ...
                    numel(value) == 1 )
                obj.Pmap_ds = value;
                obj.Pdata_valid = false;
            else
                error('??? Map has wrong dimensions. Check map size or value range.')
            end
        end
        
        function set.kf_map(obj,value)
            if isnumeric(value) && isreal(value) && ...
                    ( all( size(value) == [obj.map_size(2), obj.map_size(1)] ) || ...
                    numel(value) == 1 )
                obj.Pmap_kf = value;
            else
                error('??? Map has wrong dimensions. Check map size or value range.')
            end
        end
        
        function set.sf_map(obj,value)
            if isnumeric(value) && isreal(value) && ...
                    ( all( size(value) == [obj.map_size(2), obj.map_size(1)] ) || ...
                    numel(value) == 1 )
                obj.Pmap_sf = value;
                obj.Pdata_valid = false;
            else
                error('??? Map has wrong dimensions. Check map size or value range.')
            end
        end
        
        function set.asD_map(obj,value)
            if isnumeric(value) && isreal(value) && ...
                    ( all( size(value) == [obj.map_size(2), obj.map_size(1)] ) || ...
                    numel(value) == 1 )
                obj.Pmap_asD = value;
                obj.Pdata_valid = false;
            else
                error('??? Map has wrong dimensions. Check map size or value range.')
            end
        end
        
        function set.asA_map(obj,value)
            if isnumeric(value) && isreal(value) && ...
                    ( all( size(value) == [obj.map_size(2), obj.map_size(1)] ) || ...
                    numel(value) == 1 )
                obj.Pmap_asA = value;
                obj.Pdata_valid = false;
            else
                error('??? Map has wrong dimensions. Check map size or value range.')
            end
        end
        
        function set.esD_map(obj,value)
            if isnumeric(value) && isreal(value) && ...
                    ( all( size(value) == [obj.map_size(2), obj.map_size(1)] ) || ...
                    numel(value) == 1 )
                obj.Pmap_esD = value;
                obj.Pdata_valid = false;
            else
                error('??? Map has wrong dimensions. Check map size or value range.')
            end
        end
        
        function set.esA_map(obj,value)
            if isnumeric(value) && isreal(value) && ...
                    ( all( size(value) == [obj.map_size(2), obj.map_size(1)] ) || ...
                    numel(value) == 1 )
                obj.Pmap_esA = value;
                obj.Pdata_valid = false;
            else
                error('??? Map has wrong dimensions. Check map size or value range.')
            end
        end
        
        function set.map_valid(obj,value)
            if islogical(value)
                val = value;
            elseif all(size(value) == [1 1]) && isnumeric(value) ...
                    && isreal(value) && any(value == [0,1])
                val = logical(value);
            else
                error('??? parameter for "map_valid" is incorrect');
            end
            obj.Pmap_valid = val;
        end
        
        function set.data_valid(obj,value)
            if islogical(value)
                val = value;
            elseif all(size(value) == [1 1]) && isnumeric(value) ...
                    && isreal(value) && any(value == [0,1])
                val = logical(value);
            else
                error('??? parameter for "data_valid" is incorrect');
            end
            obj.Pdata_valid = val;
        end
                
        function set.map_extent(h_parameter_set, value)
            if ~(numel(value) == 4 && size(value, 2) == 2)
                error('Invalid value')
            end
            if ~((value(1, 2) - value(1, 1) >= 1.9*h_parameter_set.map_extension) && (value(2, 2) - value(2, 1) >= 1.9*h_parameter_set.map_extension))
                error('range should be greater than 2 times the map_extension')
            end
            h_parameter_set.Pmap_extent = value;
            
            if h_parameter_set.Psamples_per_meter > 0            
                % expand the map extent, so the map size becomes an integer
                % value without the need to round
                map_size = round(h_parameter_set.Psamples_per_meter .* (h_parameter_set.Pmap_extent(:, 2) - h_parameter_set.Pmap_extent(:, 1))+1);
                h_parameter_set.Pmap_extent(:, 2) = (map_size-1)./h_parameter_set.Psamples_per_meter + h_parameter_set.Pmap_extent(:, 1);
                h_parameter_set.Pmap_valid = false;
                h_parameter_set.Pdata_valid = false;
            end
        end
        
        function set.samples_per_meter(h_parameter_set, value)
            if ~( all(size(value) == [1 1]) && isnumeric(value) ...
                    && isreal(value) && value >= 0 )
                error('''samples_per_meter'' must be real and scalar.')
            end
            h_parameter_set.Psamples_per_meter = value;
            
            if value > 0
                % expand the map extent, so the map size becomes an integer
                % value without the need to round
                map_size = round(h_parameter_set.Psamples_per_meter .* (h_parameter_set.Pmap_extent(:, 2) - h_parameter_set.Pmap_extent(:, 1))+1);
                h_parameter_set.Pmap_extent(:, 2) = (map_size-1)./h_parameter_set.Psamples_per_meter + h_parameter_set.Pmap_extent(:, 1);
                h_parameter_set.Pmap_valid = false;
                h_parameter_set.Pdata_valid = false;
            end
        end
    end
    
    methods(Static)
        [ scenarios , config_folder , file_names ] = supported_scenarios( parse_shortnames )
    end
end

