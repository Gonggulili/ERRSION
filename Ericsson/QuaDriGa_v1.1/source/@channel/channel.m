classdef channel < handle & matlab.mixin.Copyable
%CHANNEL Class for handling channel coefficients
%
% DESCRIPTION
% Objects of this class are the output of the channel model. They are
% created by the 'channel_builder'. By default, channel
% coefficients are provided in time domain, as a list of delays and
% complex-valued amplitudes. However, this class also implements certain
% methods to post-process the channel data. Those include:     
%
%   - Transformation into frequency domain
%   - Interpolation in time domain (to change the terminal speed and
%     sampling rate) 
%   - Combining channel traces into longer segments (including birth and
%     death of clusters) 
%
% Optional PLUGIN: quadriga_channel_export
% In addition to the open-source version of QuaDRiGa, the "channel export
% plugin" provides import- and export functions for commonly used channel
% formats. This plugin adds additional methods to the "channel":
%
%    hdf5_load - Load data from stored HDF5 file
%    hdf5_save - Save data to HDF5 file
%    mport_meas_data - Convert band-limited frequency-domain data into QuaDRiGa channels
%
% QuaDRiGa Copyright (C) 2011-2015 Fraunhofer Heinrich Hertz Institute
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
    
    % Name of the 'channel' object.
    %
    %   This string is a unique identifier of the channel object. The
    %   'channel_builder' creates one channel object for each MT,
    %   each Tx and each segment. They are further grouped by scenarios
    %   (propagation environments). The string consists of four parts
    %   separated by an underscore '_'. Those are:
    %
    %       - The scenario name from 'track.scenario'
    %       - The transmitter name from 'layout.tx_name'
    %       - The receiver name from 'layout.rx_name'
    %       - The segment number
    %
    %   After 'channel.merge' has been called, the name string consists of:
    %
    %       - The transmitter name from 'layout.tx_name'
    %       - The receiver name from 'layout.rx_name'
    %
    name = 'New channel'
    
    % Version number of the QuaDRiGa release that was used to create
    % the 'channel' object.
    version = simulation_parameters.version
end

properties(Dependent)
    
    % Indicates if the path delays are identical on each MIMO link
    % (0) or if each link has different delays (1).
    individual_delays
    
    % The complex-valued channel coefficients for each path.
    %   The indices of the 4-D tensor are:
    %   [ Rx-Antenna , Tx-Antenna , Path , Snapshot ]
    coeff
    
    % The delays for each path.
    %
    %   There are two different options. If the delays are identical on
    %   the MIMO links, i.e. 'individual_delays = 0', then 'delay' is a
    %   2-D matrix with dimensions [ Path , Snapshot ].
    %
    %   If the delays are different on the MIMO links, then 'delay' is
    %   a 4-D tensor with dimensions
    %   [ Rx-Antenna , Tx-Antenna , Path , Snapshot ].
    delay
    
    
    % A flexible data structure that containing additional data
    %   This structure allows you to store additional data with the channel
    %   object. Separate variables can be identified by fieldnames 
    %   (e.g. "par.sf" or "par.snr"). Each field must contain numeric data
    %   (real or complex numbers) or strings. Nested structures are not
    %   allowed (e.g. "par.sf.val1").  
    par
    
    % The snapshot number for which the initial LSPs have been generated.
    %   Normally, this is the first snapshot. However, if the user
    %   trajectory consists of more than one segment, then
    %   'initial_position' points to the snapshot number where the
    %   current segment starts. For example: If 'initial_position' is
    %   100, then snapshots 1-99 are overlapping with the previous segment.
    initial_position
    
    % Position of each Tx in global Cartesian coordinates using units
    % of [m].
    tx_position
    
    % The receiver position global Cartesian coordinates using units of
    % [m] for each snapshot.
    rx_position
end

% Informative parameters
properties(Dependent,SetAccess=protected)
    no_rx                       % Number of receive elements (read only)
    no_tx                       % Number of transmit elements (read only)
    no_path                     % Number of paths (read only)
    no_snap                     % Number of snapshots (read only)
end

% Data storage
properties(Access=private)
    Pindividual_delays  = false;
    Pcoeff              = [];
    Pdelay              = [];
    Ppar                = [];
    Pinitial_position   = 1;
    Ptx_position        = [];
    Prx_position        = [];
end

methods
    % Constructor
    function obj = channel( Ccoeff , Cdelay , Cinitial_position , ~ )
        check = true;
        if nargin == 4
            check = false;
        end
        if nargin > 0
            if check
                obj.coeff = Ccoeff;
            else
                obj.Pcoeff = Ccoeff;
            end
        end
        if nargin > 1
            if  numel( size( Cdelay ) ) == numel( size( Ccoeff ) ) && ...
                    all( size( Cdelay ) == size( Ccoeff ) ) && ...
                    numel( size( Cdelay ) ) > 2
                obj.individual_delays = 1;
            end
            if check
                obj.delay = Cdelay;
            else
                obj.Pdelay = Cdelay;
            end
        end
        if nargin > 2
            if check
                obj.initial_position = Cinitial_position;
            else
                obj.Pinitial_position = Cinitial_position;
            end
        end
    end
    
    % Get functions
    function out = get.individual_delays(obj)
        out = obj.Pindividual_delays;
    end
    function out = get.coeff(obj)
        out = obj.Pcoeff;
    end
    function out = get.delay(obj)
        out = obj.Pdelay;
    end
    function out = get.par(obj)
        out = obj.Ppar;
    end
    function out = get.initial_position(obj)
        out = obj.Pinitial_position;
    end
    function out = get.rx_position(obj)
        out = obj.Prx_position;
    end
    function out = get.tx_position(obj)
        out = obj.Ptx_position;
    end
    function out = get.no_rx(obj)
        out = size( obj.coeff,1);
    end
    function out = get.no_tx(obj)
        out = size( obj.coeff,2);
    end
    
    function out = get.no_path(obj)
        s = size(obj.coeff);
        if numel(s) < 3
            if any(s)==0
                out = 0;
            else
                out=1;
            end
        else
            out = s(3);
        end
    end
    
    function out = get.no_snap(obj)
        s = size(obj.coeff);
        if numel(s) < 4
            if any(s)==0
                out = 0;
            else
                out=1;
            end
        else
            out = s(4);
        end
    end
    
    function v = get_version(obj)
        %GET_VERSION returns the version number
        min_ver = Inf;
        for n = 1:numel(obj)
            tmp = regexp( obj(n).version,'(?<a>[0-9]+).(?<b>[0-9]+).(?<c>[0-9]+)-(?<d>[0-9]+)' ,'tokens');
            tmp = str2double(tmp{1});
            tmp_ver = tmp(1)*1e6 + tmp(2)*1e3 + tmp(3);
            if tmp_ver < min_ver
                v = tmp;
                min_ver = tmp_ver;
            end
        end
    end
    
    % Set functions
    function set.name(obj,value)
        if ~( ischar(value) )
            error('??? "name" must be a string.')
        end
        obj.name = value;
    end
    
    function set.version(obj,value)
        if ~( ischar(value) )
            error('??? "version" must be a string.')
        elseif isempty( regexp(value,'[0-9]+.[0-9]+.[0-9]+-[0-9]+', 'once') )
            error('??? "version" must be a version-string.')
        end
        obj.version = value;
    end
    
    function set.individual_delays(obj,value)
        if ~( all(size(value) == [1 1]) && isreal(value) )
            error('??? "individual_delays" must be numeric and scalar')
        end
        value = logical( value );
        if obj.Pindividual_delays && ~value && ~isempty(obj.Pdelay)
            
            if obj.no_tx == 1 && obj.no_rx == 1
                obj.Pdelay = reshape( obj.Pdelay(1,1,:,:) , obj.no_path , obj.no_snap );
            else
                
                % Calculate power-weighted mean
                P = abs(obj.Pcoeff).^2;
                P = reshape( P, obj.no_rx*obj.no_tx , obj.no_path , obj.no_snap );
                tmp = sum( P , 1 );
                P = P ./ tmp( ones(1,obj.no_rx*obj.no_tx),:,: );
                
                D = reshape( obj.Pdelay, obj.no_rx*obj.no_tx , obj.no_path , obj.no_snap );
                D = sum( D.*P,1 );
                D( tmp == 0 ) = 0;
                D = permute( D , [2,3,1] );
                
                obj.Pdelay = D;
            end
            
        elseif ~obj.Pindividual_delays && value && ~isempty(obj.Pdelay)
            % Use the same delay for all antenna elements
            obj.Pdelay = reshape( obj.Pdelay,1,1,obj.no_path,obj.no_snap );
            obj.Pdelay = obj.Pdelay( ones(1,obj.no_rx) , ones(1,obj.no_tx) , : , :  );
        end
        obj.Pindividual_delays = value;
    end
    
    function set.coeff(obj,value)
        if ~isnumeric(value)
            error('??? "coeff" must be numeric')
            
        elseif isempty( value ) % Remove coefficients
            obj.Pcoeff = [];
            obj.Pdelay = [];
            obj.Prx_position = [];
            obj.Pinitial_position = 0;
            
        else
            so = size(obj.Pcoeff);
            sn = size(value);
            if (numel(so) == 2 && numel(sn) == 2) ||...
                    (numel(so) >= 3 && numel(sn) >= 3 && so(3) == sn(3))
                obj.Pcoeff = value;
            else
                obj.Pcoeff = value;
                if obj.Pindividual_delays
                    obj.Pdelay = zeros( obj.no_rx , obj.no_tx , obj.no_path, obj.no_snap );
                else
                    obj.Pdelay = zeros( obj.no_path , obj.no_snap );
                end
            end
        end
    end
    
    function set.delay(obj,value)
        [a, b, c, d] = size(value);
        if ~isnumeric(value)
            error('??? "delay" must be numeric')
        elseif obj.individual_delays && ...
                ~( all( [a,b,c,d] == [ obj.no_rx , obj.no_tx , obj.no_path , obj.no_snap ] ) || ...
                any( [a,b,c,d]==0 ) )
            error('??? "delay" must match the number of antennas, paths and snapshots')
        elseif ~obj.individual_delays && ...
                ~all( [a,b] == [ obj.no_path , obj.no_snap ] )
            error('??? "delay" must match the number of paths and snapshots')
        end
        obj.Pdelay = value;
    end
    
    function set.par(obj,value)
        if ~isempty( value )
            if ~isstruct( value )
                error('??? "par" must be a structure')
            else
                names = fieldnames( value );
                for i_names = 1:numel( names )
                    if isstruct( value.( names{i_names} ) )
                        error('??? "par" cannot contain nested structures');
                    elseif iscell( value.( names{i_names} ) )
                        error('??? "par" cannot contain cell arrays');
                    end
                end
                obj.Ppar = value;
            end
        else
            obj.Ppar = [];
        end
    end
    
    function set.initial_position(obj,value)
        if ~( all(size(value) == [1 1]) && isnumeric(value) ...
                && isreal(value) && mod(value,1)==0 && value >= 0 )
            error('??? "initial_position" must be integer and >= 0')
        elseif value > obj.no_snap
            error('??? "initial_position" must not exceed the number of snapshots')
        end
        obj.Pinitial_position = value;
    end
    
    function set.rx_position(obj,value)
        if isempty( value)
            obj.Prx_position = [];
        else
            if ~( isnumeric(value) && isreal(value) && size(value,2) > 0 )
                error('QuaDRiGa:Channel:wrongInputValue','??? "rx_position" must consist of real numbers')
            elseif ~all( size(value,1) == 3 )
                error('QuaDRiGa:Channel:wrongInputValue','??? "rx_position" must have three rows.')
            elseif size( value,2 ) ~= obj.no_snap
                error('QuaDRiGa:Channel:wrongInputValue','??? "rx_position" must match the number of snapshots.')
            end
            obj.Prx_position = value;
        end
    end
    
    function set.tx_position(obj,value)
        if isempty( value)
            obj.Ptx_position = [];
        else
            if ~( isnumeric(value) && isreal(value) )
                error('QuaDRiGa:Channel:wrongInputValue','??? "tx_position" must consist of real numbers')
            elseif ~all( size(value,1) == 3 && size(value,2) == 1 )
                error('QuaDRiGa:Channel:wrongInputValue','??? "tx_position" must have 3 rows')
            end
            obj.Ptx_position = value;
        end
    end
end

methods(Static)
    [ obj,dims ] = hdf5_load( varargin )
    [ h_channel, snr ] = import_meas_data( Y, B, L_max, usage, noise_limit, delay_limit, pilot_grid, verbose, show_pdp )
end

end
