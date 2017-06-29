classdef array < handle & matlab.mixin.Copyable
%ARRAY Antenna array class
%
% DESCRIPTION
% This class combines all functions to create and edit antenna arrays. An antenna
% array is a set of single antenna elements, each having a specific beam pattern,
% that can be combined in any geometric arrangement. A set of synthetic arrays
% that allow simulations without providing your own antenna patterns is
% provided (see generate method for more details).
%
% REFERENCE
% The main functionality was taken from the Winner channel model. "Kyösti, P.;
% Meinilä, J.; Hentilä, L. & others; {IST-4-027756 WINNER II D1.1.2 v.1.1}:
% WINNER II Channel Models; 2007". New functionality has been added to provide
% geometric polarization calculations and antenna pattern rotations.
%
% EXAMPLE
% This example creates an array of crossed dipoles.
%
%    a = array;                           % Create new array object
%    a.generate('dipole');                % Generate a synthetic dipole pattern
%    a.copy_element(1,2);                 % Duplicate the dipole
%    a.rotate_pattern(90,'y',2);          % Rotate the second element by 90°
%    a.visualize;                         % Show the output
%
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
        name = 'New array';                     % Name of the antenna array
        
        % Method for interpolating the beam patterns
        %   The default is linear interpolation. Optional are:
        %       nearest - Nearest neighbor interpolation (QuaDRiGa optimized)
        %       linear - Linear interpolation (QuaDRiGa optimized, Default)
        %       spline - Cubic spline interpolation (MATLAB internal function)
        %       nearest_int - Nearest neighbor interpolation (MATLAB internal function)
        %       linear_int - Linear interpolation (MATLAB internal function)
        %
        %   Note: MATLAB internal routines slow down the simulations significantly.
        interpolation_method = 'linear';        
        
        % The polarization basis of the pattern
        %   The following bases are currently supported:
        %
        %       cartesian           - Ludwig 1
        %       az-el               - Ludwig 2 - Azimuth over Elevation
        %       el-az               - Ludwig 2 - Elevation over Azimuth
        %       polar-spheric       - Ludwig 2 - Polar-Spheric [DEFAULT]
        %
        %   You can specify the polarization basis of the pattern by setting
        %   the appropriate string. By default, QuaDRiGa requires a
        %   polar-spheric basis. If a different basis is specified, an
        %   appropriate transformation will be carried out.
        polarization_basis = 'polar-spheric';
    end
    
    properties(Dependent)
        % Numeric precision
        %   Single precision arrays need only half the storage space and
        %   computations are usually faster.
        %   Default: Double precision
        precision
        
        % Number of antenna elements in the array
        %    Increasing the number of elements creates new elements which
        %    are initialized as copies of the first element. Decreasing the
        %    number of elements deletes the last elements from the array.
        no_elements
        
        % Elevation angles in [rad] were samples of the field patterns are provided
        %   The field patterns are given in spherical coordinates. This
        %   variable provides the elevation sampling angles in radians
        %   ranging from -pi/2 (downwards) to pi/2 (upwards).  
        elevation_grid

        % Azimuth angles in [rad] were samples of the field patterns are provided
        %   The field patterns are given in spherical coordinates. This
        %   variable provides the azimuth sampling angles in radians
        %   ranging from -pi to pi.  
        azimuth_grid
        
        % Position of the antenna elements in local Cartesian coordinates
        % (units of [m])
        element_position
        
        % The first component of the antenna pattern. If the polar-spheric
        % polarization basis is used, this variable contains the vertical
        % (or theta) component of the electric field given in spherical
        % coordinates.   
        %   This variable is a tensor with dimensions [ elevation, azimuth,
        %   element ] describing the vertical (or theta) component of the
        %   far field of each antenna element in the array.  
        Fa
        
        % The second component of the antenna pattern. If the polar-spheric
        % polarization basis is used, this variable contains the horizontal
        % (or phi) component of the electric field given in spherical
        % coordinates.   
        %   This variable is a tensor with dimensions [ elevation, azimuth,
        %   element ] describing the horizontal (or phi) component of the
        %   far field of each antenna element in the array.  
        Fb
        
        % The third component of the antenna pattern. Currently, it is only
        % used when the antenna pattern is using a Cartesian polarization
        % basis.  
        Fc
        
        % Coupling matrix between elements
        %   This matrix describes a pre- or post-processing of the signals
        %   that are fed to the antenna elements. For example, in order to
        %   transmit a LHCP signal, two antenna elements are needed.
        %   They are then coupled by a matrix
        %
        %   	1/sqrt(2) * [1;j]
        %
        %   The rows in the matrix correspond to the antenna elements, the
        %   columns to the signal ports. In this example, the antenna has
        %   one port, i.e. it is fed with one input signal. This signal is
        %   then split into two and fed to the two antenna elements where
        %   the second element radiates the signal with 90 degree phase
        %   shift.
        %   In a similar fashion, it is possible to create fixed
        %   beamforming antennas and include crosstalk between antenna
        %   elements. By default, coupling is set to an identity matrix
        %   which indicates perfect isolation between the antenna elements.           
        coupling
    end

    properties(Dependent,SetAccess=private)
        % Indicates if the array is compressed
        %    It is possible to compress the antenna array in memory to save
        %    storage space and relay the memory requirements for large
        %    arrays. This property indicates it the array is compressed.
        iscompressed
    end
    
    % The legacy handles for the field pattern
    properties(Dependent,Hidden)
        field_pattern_vertical
        field_pattern_horizontal
    end
    
    properties(SetAccess=private)
        no_az = 5;                              % Number of azimuth values
        no_el = 3;                              % Number of elevation values
    end

    properties(Access=private)
        Pprecision                  = 'double';
        Pno_elements                = 1;
        Pelevation_grid             = [ -1.570796326794897,0,1.570796326794897];
        Pazimuth_grid               = [ -3.141592653589793,-1.570796326794897,0,...
            1.570796326794897,3.141592653589793];
        Pelement_position           = [0;0;0];
        PFa                         = ones(3,5);
        PFb                         = zeros(3,5);
        PFc                         = 0;
        Pind                        = [];
        Pcoupling                   = 1;
    end

    methods
        % The constructor
        function h_array = array( array_type, varargin )
            if nargin > 0
                h_array.generate( array_type , 0 , varargin{:} );
            else
                h_array.generate( 'omni' );
            end
        end
        
        % Get functions
        function out = get.precision(obj)
            out = obj.Pprecision;
        end
        function out = get.iscompressed(obj)
            out = ~isempty(obj.Pind);
        end
        function out = get.no_elements(obj)
            out = obj.Pno_elements;
        end
        function out = get.elevation_grid(obj)
            out = obj.Pelevation_grid;
        end
        function out = get.azimuth_grid(obj)
            out = obj.Pazimuth_grid;
        end
        function out = get.element_position(obj)
            out = obj.Pelement_position;
        end
        function out = get.Fa(obj)
            if obj.iscompressed
                out = obj.PFa( :,:, obj.Pind(1,:) );
            else
                out = obj.PFa;
            end
        end
        function out = get.Fb(obj)
            if obj.iscompressed
                out = obj.PFb( :,:, obj.Pind(2,:) );
            else
                out = obj.PFb;
            end
        end
        function out = get.Fc(obj)
            out = obj.PFc;
        end
        function out = get.coupling(obj)
            out = obj.Pcoupling;
        end
        function out = get.field_pattern_vertical(obj)
            out = obj.PFa;
        end
        function out = get.field_pattern_horizontal(obj)
            out = obj.PFb;
        end
        
        % Set functions
        function set.name(obj,value)
            if ~( ischar(value) )
                error('QuaDRiGa:Array:wrongInputValue','??? "name" must be a string.')
            end
            obj.name = value;
        end
        
        function set.interpolation_method(obj,value)
            supported_types = {'nearest','linear','spline','linear_int','nearest_int'};
            if ~( ischar(value) && any( strcmpi(value,supported_types)) )
                str = 'Interpolation method not supported; supported types are: ';
                no = numel(supported_types);
                for n = 1:no
                    str = [str,supported_types{n}];
                    if n<no
                        str = [str,', '];
                    end
                end
                error('QuaDRiGa:Array:wrongInterpolationType',str);
            end
            obj.interpolation_method = value;
        end
        
        function set.polarization_basis(obj,value)
            if ~array.supported_pol_basis( value )
                str = 'Polarization basis not supported.';
                error('QuaDRiGa:Array:wrongPolarizationBasis',str);
            end
            obj.polarization_basis = value;
        end
        
        function set.precision(obj,value)
            supported_types = {'double','single'};
            if ~( ischar(value) && any( strcmpi(value,supported_types)) )
                str = 'Precision not supported; supported types are: ';
                no = numel(supported_types);
                for n = 1:no
                    str = [str,supported_types{n}];
                    if n<no
                        str = [str,', '];
                    end
                end
                error('QuaDRiGa:Array:wrongPrecision',str);
            end
            
            switch value
                case 'double'
                    obj.Pelevation_grid = double( obj.Pelevation_grid );
                    obj.Pazimuth_grid = double( obj.Pazimuth_grid );
                    obj.Pelement_position = double( obj.Pelement_position );
                    obj.PFa = double( obj.PFa );
                    obj.PFb = double( obj.PFb );
                    obj.PFc = double( obj.PFc );
                    obj.Pcoupling  = double( obj.Pcoupling );
                    
                case 'single'
                    obj.Pelevation_grid = single( obj.Pelevation_grid );
                    obj.Pazimuth_grid = single( obj.Pazimuth_grid );
                    obj.Pelement_position = single( obj.Pelement_position );
                    obj.PFa = single( obj.PFa );
                    obj.PFb = single( obj.PFb );
                    obj.PFc = single( obj.PFc );
                    obj.Pcoupling  = single( obj.Pcoupling );
                    
            end
            obj.Pprecision = value;
        end
        
        function set.no_elements(obj,value)
            if ~( all(size(value) == [1 1]) && isnumeric(value) ...
                    && isreal(value) && mod(value,1)==0 && value > 0 )
                error('QuaDRiGa:Array:wrongInputValue','??? "no_elements" must be integer and > 0')
            end
            
            % Requires uncompressed arrays
            obj.uncompress;
            
            if obj.no_elements > value
                obj.Pelement_position = obj.Pelement_position(:,1:value);
                obj.PFa = obj.PFa(:,:,1:value);
                obj.PFb = obj.PFb(:,:,1:value);
                
                if numel( obj.PFc ) ~= 1
                    obj.PFc = obj.PFc(:,:,1:value);
                end
                                
                ne = obj.no_elements-value;
                nc = size(obj.Pcoupling);
                obj.Pcoupling = obj.Pcoupling( 1:nc(1)-ne , 1:max(nc(2)-ne,1) );
                
            elseif obj.no_elements < value
                ne = value-obj.no_elements;
                
                obj.Pelement_position = [ obj.Pelement_position,...
                    obj.Pelement_position(:,ones( 1,ne )) ];
                
                obj.PFa = cat( 3, obj.PFa ,...
                    obj.PFa(:,:,ones( 1,ne )));
                
                obj.PFb = cat( 3, obj.PFb ,...
                    obj.PFb(:,:,ones( 1,ne )));
                
                if numel( obj.PFc ) ~= 1
                    obj.PFc = cat( 3, obj.PFc ,...
                        obj.PFc(:,:,ones( 1,ne )));
                end
                   
                nc = size(obj.Pcoupling);
                C = zeros( nc(1)+ne , nc(2)+ne);
                for n = 1:ne
                   C( nc(1)+n,nc(2)+n ) = 1; 
                end
                C( 1:nc(1) , 1:nc(2) ) = obj.Pcoupling;
                obj.Pcoupling = C;
            end
            
            obj.Pno_elements = value;
        end
        
        function set.elevation_grid(obj,value)
            if ~( any( size(value) == 1 ) && isnumeric(value) && isreal(value) &&...
                    max(value)<=pi/2+1e-7 && min(value)>=-pi/2-1e-7 )
                error('QuaDRiGa:Array:wrongInputValue','??? "elevation_grid" must be a vector containing values between -pi/2 and pi/2')
            end
            
            % Perform typecast
            switch obj.Pprecision
                case 'double'
                    value = double( value );
                case 'single'
                    value = single( value );
            end
            
            % Requires uncompressed arrays
            obj.uncompress;
            
            val_old = size(obj.Fa,1);
            val_new = numel(value);
            if val_old > val_new
                obj.PFa = obj.PFa( 1:val_new ,: , : );
                obj.PFb = obj.PFb( 1:val_new ,: , : );
                if numel( obj.PFc ) ~= 1
                    obj.PFc = obj.PFc( 1:val_new ,: , : );
                end
                
            elseif val_old < val_new
                a = size( obj.Fa );
                if numel(a) == 2
                    a(3) = 1;
                end
                b = val_new-val_old;
                
                obj.PFb = cat( 1, obj.PFb ,...
                    ones( b ,a(2),a(3))  );
                
                obj.PFa = cat( 1, obj.PFa ,...
                    ones( b ,a(2),a(3))  );
                
                if numel( obj.PFc ) ~= 1
                    obj.PFc = cat( 1, obj.PFc ,...
                        ones( b ,a(2),a(3))  );
                end
            end
            
            if size(value,1) ~= 1
                obj.Pelevation_grid = value';
            else
                obj.Pelevation_grid = value;
            end
            obj.no_el = val_new;
        end
        
        function set.azimuth_grid(obj,value)
            if ~( any( size(value) == 1 ) && isnumeric(value) && isreal(value) &&...
                     max(value)<=pi+1e-7 && min(value)>=-pi-1e-7 )
                error('QuaDRiGa:Array:wrongInputValue','??? "azimuth_grid" must be a vector containing values between -pi and pi')
            end
            
            % Perform typecast
            switch obj.Pprecision
                case 'double'
                    value = double( value );
                case 'single'
                    value = single( value );
            end
            
            % Requires uncompressed arrays
            obj.uncompress;
            
            val_old = size(obj.Fa,2);
            val_new = numel(value);
            if val_old > val_new
                obj.PFb = obj.PFb( : , 1:val_new , : );
                obj.PFa = obj.PFa( : , 1:val_new  , : );
                if numel( obj.PFc ) ~= 1
                    obj.PFc = obj.PFc( : , 1:val_new  , : );
                end
                
            elseif val_old < val_new
                a = size( obj.Fa );
                if numel(a) == 2
                    a(3) = 1;
                end
                b = val_new-val_old;
                
                obj.PFb = cat( 2, obj.PFb ,...
                    ones( a(1) , b , a(3))  );
                
                obj.PFa = cat( 2, obj.PFa ,...
                    ones( a(1) , b , a(3))  );
                
                if numel( obj.PFc ) ~= 1
                    obj.PFc = cat( 1, obj.PFc ,...
                        ones( a(1) , b , a(3))  );
                end
            end
            
            if size(value,1) ~= 1
                obj.Pazimuth_grid = value';
            else
                obj.Pazimuth_grid = value;
            end
            obj.no_az = val_new;
        end
        
        function set.element_position(obj,value)
            if ~( isnumeric(value) && isreal(value) )
                error('QuaDRiGa:Array:wrongInputValue','??? "element_position" must consist of real numbers')
            elseif ~all( size(value,1) == 3 )
                error('QuaDRiGa:Array:wrongInputValue','??? "element_position" must have 3 rows')
            end
            if size(value,2) ~= obj.no_elements
                obj.no_elements = size(value,2);
            end
            
            % Perform typecast
            switch obj.Pprecision
                case 'double'
                    value = double( value );
                case 'single'
                    value = single( value );
            end
            
            obj.Pelement_position = value;
        end
        
        function set.Fa(obj,value)
            a = numel( obj.Pelevation_grid );
            b = numel( obj.Pazimuth_grid );
            
            if obj.Pno_elements == 1
                dims = [ a , b ];
            else
                dims = [ a , b , obj.Pno_elements];
            end
            
            if ~( isnumeric(value) )
                error('QuaDRiGa:Array:wrongInputValue','??? "Fa" must be numeric.')
            elseif ~( numel(size(value)) == numel(dims) && all( size(value) == dims ) )
                error('QuaDRiGa:Array:wrongInputValue',['??? "Fa" must be of size [',num2str(a),'x',num2str(b),...
                    'x',num2str(obj.Pno_elements),'].'])
            end
            
            % Requires uncompressed arrays
            obj.uncompress;
            
            % Perform typecast
            switch obj.Pprecision
                case 'double'
                    value = double( value );
                case 'single'
                    value = single( value );
            end
            
            obj.PFa = value;
        end
        
        function set.Fb(obj,value)
            a = numel( obj.Pelevation_grid );
            b = numel( obj.Pazimuth_grid );
            
            if obj.no_elements == 1
                dims = [ a , b ];
            else
                dims = [ a , b , obj.Pno_elements];
            end
            
            if ~( isnumeric(value) )
                error('QuaDRiGa:Array:wrongInputValue','??? "Fb" must be numeric.')
            elseif ~( numel( size(value) ) == numel( dims ) && all( size(value) == dims ) )
                error('QuaDRiGa:Array:wrongInputValue',['??? "Fb" must be of size [',num2str(a),'x',num2str(b),...
                    'x',num2str(obj.Pno_elements),'].'])
            end
            
            % Requires uncompressed arrays
            obj.uncompress;
            
            % Perform typecast
            switch obj.Pprecision
                case 'double'
                    value = double( value );
                case 'single'
                    value = single( value );
            end
            
            obj.PFb = value;
        end
        
        function set.Fc(obj,value)
            if numel( value ) == 1 && value == 0
                obj.PFc = 0;
            else
                a = numel( obj.Pelevation_grid );
                b = numel( obj.Pazimuth_grid );
                
                if obj.no_elements == 1
                    dims = [ a , b ];
                else
                    dims = [ a , b , obj.Pno_elements];
                end
                
                if ~( isnumeric(value) )
                    error('QuaDRiGa:Array:wrongInputValue','??? "Fc" must be numeric.')
                elseif ~( numel( size(value) ) == numel( dims ) && all( size(value) == dims ) )
                    error('QuaDRiGa:Array:wrongInputValue',['??? "Fc" must be of size [',num2str(a),'x',num2str(b),...
                        'x',num2str(obj.Pno_elements),'].'])
                end
                
                % Requires uncompressed arrays
                obj.uncompress;
                
                % Perform typecast
                switch obj.Pprecision
                    case 'double'
                        value = double( value );
                    case 'single'
                        value = single( value );
                end
                
                obj.PFc = value;
            end
        end
       
        function set.coupling(obj,value)
            if ~( isnumeric(value) && size(value,1) == obj.Pno_elements && ...
                    size(value,2) >= 1 )
                error('QuaDRiGa:Array:wrongInputValue','??? "coupling" must be a matrix with rows equal to elements and columns equal to ports')
            end
            
            % Perform typecast
            switch obj.Pprecision
                case 'double'
                    value = double( value );
                case 'single'
                    value = single( value );
            end
            
            obj.Pcoupling = value;
        end
        
        % Handles for the old field-pattern parameters
        function set.field_pattern_vertical(obj,value)
            obj.Fa = value;
        end
        function set.field_pattern_horizontal(obj,value)
            obj.Fb = value;
        end
    end

    methods(Static)
        function types = supported_types
            types =  {'omni', 'dipole', 'short-dipole', 'half-wave-dipole',...
                'custom', '3gpp-macro','3gpp-3d',...
                'parametric','rhcp-dipole', 'lhcp-dipole', 'lhcp-rhcp-dipole', ...
                'xpol', 'ula2', 'ula4', 'ula8', 'patch'};
        end
        h_array = import_pattern( fVi, fHi )
        out = supported_pol_basis( value )
    end
end
