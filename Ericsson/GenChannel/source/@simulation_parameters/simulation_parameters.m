classdef simulation_parameters < handle & matlab.mixin.Copyable
%SIMULATION_PARAMETERS General configuration settings
%
% DESCRIPTION
% This class controls the simulation options and calculates constants for other
% classes. Currently, the following options can be set:
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
    
    properties(Dependent)
        
        % The number of samples per half-wave length
        % 	Sampling density describes the number of samples per half-wave
        % 	length. To fulfill the sampling theorem, the minimum sample
        % 	density must be 2. For  smaller values, interpolation of the
        % 	channel for variable speed is not  possible. On the other hand,
        % 	high values significantly increase the computing time
        % 	significantly. A good value is around 4. 
        sample_density                              
        samples_per_meter                           % Samples per one meter
        
        % Precision of the drifting functionality
        %   drifting_precision = 0
        %   This method applies rotating phasors to each path which
        %   emulates time varying Doppler characteristics. However, the
        %   large-scale parameters (departure and arrival angles, shadow
        %   fading, delays, etc.) are not updated in this case. This mode
        %   requires the least computing resources and may be preferred
        %   when only short linear tracks (up to several cm) are considered
        %   and the distance between transmitter and receiver is large. The
        %   phases at the antenna arrays are calculated by a planar wave
        %   approximation.
        %
        %   drifting_precision = 1 (default)
        %   When drifting is enabled, all arrival angles, the LOS departure
        %   angle, delays, and phases are updated for each snapshot using a
        %   single-bounce model. This requires significantly more computing
        %   resources but also increases the accuracy of the results.
        %   Drifting is required when non-linear tracks are generated or
        %   the distance between transmitter and receiver is small (below
        %   20 m). The phases at the antenna arrays are calculated by a
        %   planar wave approximation.
        %
        %   drifting_precision = 2
        %   The arrival angles, the LOS departure angle, delays, and phases
        %   are updated for each snapshot and for each antenna element at
        %   the receiver (spherical wave assumption). The phases at the
        %   transmitter are calculated by a planar wave approximation. This
        %   increases the accuracy for multi-element antenna arrays at the
        %   receiver. However, the computational complexity increases as
        %   well.
        %
        %   drifting_precision = 3 [EXPERIMENTAL]
        %   This option also calculates the shadow fading, path loss and
        %   K-factor for each antenna element at the receiver separately.
        %   This feature tends to predict higher MIMO capacities since is
        %   also increases the randomness of the power for different MIMO
        %   elements.    
        %
        %   drifting_precision = 4
        %   This option uses spherical waves at both ends, the transmitter
        %   and the receiver. This method assumes a single-bounce model and
        %   no mapping of departure and arrival angles is done. Hence,
        %   departure angular spreads are effectively ignored and results
        %   might be erroneous.
        %
        %   drifting_precision = 5 [EXPERIMENTAL]
        %   This option uses spherical waves at both ends, the transmitter
        %   and the receiver. This method uses a multi-bounce model where
        %   the departure and arrival angels are matched such that the
        %   angular spreads stay consistent.
        drifting_precision
        
    end
    
    properties
        % Select the polarization rotation method
        %   use_polarization_rotation = 0
        %   Uses the polarization method from WINNER. No polarization
        %   rotation is calculated.
        %
        %   use_polarization_rotation = 1
        %   Uses the new polarization rotation method where the XPR is
        %   modeled by a rotation matrix. No change of circular
        %   polarization is assumed. 
        %
        %   use_polarization_rotation = 2 (default)
        %   Uses the polarization rotation with an additional phase offset
        %   between the H and V component of the NLOS paths. The offset
        %   angle is calculated to match the XPR for circular polarization.
        %
        %   use_polarization_rotation = 3
        %   Uses polarization rotation for the geometric polarization but
        %   models the NLOS polarization change as in WINNER.
        use_polarization_rotation   = 2;            
        
        % Returns absolute delays in channel impulse response 
        %   By default, delays are calculated such that the LOS delay is
        %   normalized to 0. By setting use_absolute_delays to 1 or true,
        %   the absolute path delays are included in channel.delays at the
        %   output of the model.
        use_absolute_delays         = false;
        
        % Generates an additional path for the ground reflection
        %   If this variable is set to true, the channel model will create an additional
        %   path for the ground reflections. The basic assumption is that the ground
        %   reflection is not delay-resolvable and thus can be considered part of the LOS
        %   path. Hence, the LOS path gets split into two paths, the direct path and the
        %   ground reflection.
        use_ground_reflection       = false;
        
        % Initializes each path with a random initial phase
        %   By default, each path is initialized with a random phase.
        %   Setting "use_random_initial_phase" to zeros disables this
        %   function. In this case, each path gets initialized with a
        %   zero-phase.
        use_random_initial_phase    = true;
        
        % Selects the angular mapping method
        %   use_angular_mapping = 1
        %   Maps the path powers to arrival angles by a wrapped Gaussian
        %   distribution. This method is adopted from the WINNER model.
        %   However, the generated angles show high correlations if the
        %   K-Factor is larger than 0 dB.
        %
        %   use_angular_mapping = 2 [Default]
        %   This method generates random angles for the paths. The angular
        %   spread is maintained by a scaling operation. The output angles
        %   have a more natural distribution. However, there is an upper
        %   limit for the angular spread of roughly 100 degree in NLOS
        %   conditions.
        use_angular_mapping         = 2;
        
        % Selects the parameter map generation algorithm
        %   use_map_algorithm = 1
        %   Uses the algorithm from the WINNER model. 
        %   Bakowski, K. & Wesolowski, K.; Change the Channel; IEEE Veh.
        %   Technol. Mag., 2011, 6, 82-91 
        %
        %   use_map_algorithm = 2 [Default]
        %   Uses a modified version of the WINNER algorithm that also
        %   filters the diagonal directions.
        use_map_algorithm           = 2;
        
        show_progress_bars          = true;         % Show a progress bar on the MATLAB prompt
        center_frequency            = 2.6e9;        % Center frequency in [Hz]
        
        % Resolution of the LSP maps in [samples/m]
        % Setting a value of 0 automatically chooses the optimal map
        % resolution depending on the values in the parameter-tables.
        map_resolution              = 0;
    end
    
    properties(Constant)
        version = '1.4.8-571';    % Version number of the current QuaDRiGa release (constant)
    end
    
    properties(Constant)
        speed_of_light = 299792458;                 % Speed of light (constant)
    end
    
    properties(Dependent,SetAccess=protected)
        wavelength                                  % Carrier wavelength in [m] (read only)
    end
    
    properties(Access=private)
        Psample_density             = 2.5;
        Pdrifting_precision         = 1;
    end
    
    methods
        
        % Get functions
        function out = get.sample_density(obj)
            out = obj.Psample_density;
        end
        function out = get.samples_per_meter(obj)
            out = 2*obj.center_frequency*obj.Psample_density / obj.speed_of_light;
        end
        function out = get.wavelength(obj)
            out = obj.speed_of_light / obj.center_frequency;
        end
        function out = get.drifting_precision(obj)
            out = obj.Pdrifting_precision;
        end
        
        % Set functions
        function set.sample_density(obj,value)
            if ~( all(size(value) == [1 1]) && isnumeric(value) && isreal(value) && value > 0 )
                error('??? Invalid sample density. The value must be real and > 0.')
            end
            obj.Psample_density = value;
        end
        
        function set.samples_per_meter(obj,value)
            if ~( all(size(value) == [1 1]) && isnumeric(value) && isreal(value) && value > 0 )
                error('??? Invalid samples_per_meter. The value must be real and > 0.')
            end
            obj.Psample_density = value*obj.wavelength/2;
        end
        
        function set.center_frequency(obj,value)
            if ~( all(size(value) == [1 1]) && isnumeric(value) && isreal(value) && value >= 0 )
                error('??? Invalid center frequency. The value must be real and > 0.')
            end
            obj.center_frequency = value;
        end
        
        function set.use_absolute_delays(obj,value)
            if ~( all(size(value) == [1 1]) ...
                    && (isnumeric(value) || islogical(value)) ...
                    && any( value == [0 1] ) )
                error('??? "use_absolute_delays" must be 0 or 1')
            end
            obj.use_absolute_delays = logical( value );
        end
        
        function set.use_ground_reflection(obj,value)
            if ~( all(size(value) == [1 1]) ...
                    && (isnumeric(value) || islogical(value)) ...
                    && any( value == [0 1] ) )
                error('??? "use_ground_reflection" must be 0 or 1')
            end
            obj.use_ground_reflection = logical( value );
        end
        
        function set.drifting_precision(obj,value)
            if ~( all(size(value) == [1 1]) ...
                    && isnumeric(value) ...
                    && any( value == 0:5 ) )
                error(['??? "drifting_precision" must be in between 0 and 5.',...
                    ' Type ''help simulation_parameters.drifting_precision'' for more infos.'])
            end
            obj.Pdrifting_precision = value;
        end
        
        function set.use_map_algorithm(obj,value)
            if ~( all(size(value) == [1 1]) ...
                    && isnumeric(value) ...
                    && any( value == 1:2 ) )
                error('??? "use_map_algorithm" must be 1 or 2.');
            end
            obj.use_map_algorithm = value;
        end
        
        function set.use_polarization_rotation(obj,value)
            if ~( all(size(value) == [1 1]) ...
                    && (isnumeric(value) || islogical(value)) ...
                    && any( value == 0:3 ) )
                error('??? "use_polarization_rotation" must be 0, 1, 2 or 3')
            end
            obj.use_polarization_rotation = value;
        end
        
        function set.map_resolution(obj,value)
            if ~( all(size(value) == [1 1]) && isnumeric(value) ...
                    && isreal(value) && value >= 0 )
                error('??? "map_resolution" must be real and scalar.')
            end
            obj.map_resolution = value;
        end
        
        function set.show_progress_bars(obj,value)
            if ~( all(size(value) == [1 1]) ...
                    && (isnumeric(value) || islogical(value)) ...
                    && any( value == [0 1] ) )
                error('??? "use_subpath_output" must be 0 or 1')
            end
            obj.show_progress_bars = logical( value );
        end
    end
end

