classdef channel_builder < handle
%CHANNEL_BUILDER Class for generating the channel coefficients
%
% DESCRIPTION
% This class implements all functions that are needed to generate the
% channel coefficients. It thus implements the core components of the
% channel model. The class holds all the input variables as properties.
% It's main function 'get_channels' then generates the coefficients. The
% procedure is summarized as follows:
%
% The channel builder first generates a set of random clusters around each
% receiver. This is done by drawing random variables for the delay, the
% power and the departure and arrival angles for each cluster. Each cluster
% thus represents the origin of a reflected (and scattered) signal. The
% clusters are then represented as taps in the final CIR. The random
% variables fit the distributions and correlations defined by the
% parameter_set object.
%
% Next, antenna dependent parameters are extracted for each user. Those
% depend on the position of the terminal, its orientation and the equipped
% antennas. The polarization rotation of the NLOS taps is modeled by a
% random variable which fits to the  distribution defined by the
% parameter_set. The LOS polarization is calculated from the geometric
% orientation of the antennas. A core function here is the interpolation of
% the antenna patterns which results in a specific H and V value for each
% subpath.
%
% The core function then generates the coefficients themselves. This is
% done for each antenna element and for each snapshot separately and
% also includes the Doppler shift of each subpath. Finally, the k-factor
% and the shadow fading are applied and a all the data is returned as an
% channel object.
%
% REFERENCE
% Some functions were taken from the Winner channel model. "Kyösti, P.;
% Meinilä, J.; Hentilä, L. & others; {IST-4-027756 WINNER II D1.1.2 v.1.1}:
% WINNER II Channel Models; 2007". New functions support the geometric
% generation of polarized channels and non-linear tracks.
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
        name = 'channel';    	% Name of the 'channel_builder' object
        par                    	% The 'parameter_set' object for this channel builder
        
        % The number of clusters. 
        NumClusters
        
        % The initial delays for each path in [s]. Rows correspond to the
        % MTs, columns to the paths. 
        taus
        
        % The normalized initial power (squared average amplitude) for each
        % path. Rows correspond to the MT, columns to the paths. The sum
        % over all columns must be 1.   
        pow
        
        AoD                    	% The initial azimuth of departure angles for each path in [rad].
        AoA                    	% The initial azimuth of arrival angles for each path in [rad].
        EoD                   	% The initial elevation of departure angles for each path in [rad].
        EoA                    	% The initial elevation of departure angles for each path in [rad].
        
        % The initial cross polarization power ratio in [dB] for each
        % sub-path. The dimensions correspond to the MT, the path number,
        % and the sub-path number.   
        xpr
        
        pin                   	% The initial phases in [rad] for each sub-path.
        
        % The phase offset angle for the circular XPR in [rad]. The
        % dimensions correspond to the MT, the path number, and the
        % sub-path number.   
        kappa
        
        % Random phasors for the WINNER polarization coupling method.
        % The dimensions correspond to polarization matrix index
        % '[ 1 3 ; 2 4 ]}', the subpath number and the MT.
        random_pol            	% Random Polarization matrix
        
        % A random index list for the mutual coupling of subpaths at the Tx
        % and Rx. The dimensions correspond to the subpath index (1-20),
        % the angle (AoD, AoA, EoD, EoA), the path number and the MT.
        subpath_coupling
    end

    properties(Access = private)
        pow_wo_kf
    end
    
    methods
        % Constructor
        function h_cb = channel_builder( h_parset )
            if numel( h_parset ) > 1
                error('"h_parset" must be scalar')
            end
            h_cb.name        = h_parset.name;
            h_cb.par         = h_parset;
        end
        
        % Set-Functions
        function set.name(obj, value)
            if ~(ischar(value))
                error('??? "name" must be a string.')
            end
            obj.name = value;
        end
        
        function set.par(obj, value)
            if ~(isa(value, 'parameter_set'))
                error('??? "par" must be of class "parameter_set".')
            elseif ~all(size(value) == [1, 1])
                error('??? "par" must be scalar.')
            end
            
            if value.no_positions > 1 && numel( value.rx_array ) == 1
                value.rx_array = value.rx_array(ones(1,value.no_positions));
            end
            if size(value.rx_array, 1) ~= 1
                value.rx_array = value.rx_array';
            end
            
            if value.no_positions > 1 && numel( value.rx_track ) == 1
                value.rx_track = value.rx_track(ones(1,value.no_positions));
            end
            if size(value.rx_track, 1) ~= 1
                value.rx_track = value.rx_track';
            end
            
            for n = 1:value.no_positions
                if ~(value.rx_track(n).no_segments <= 2)
                    error('??? Each "track" must have at most 2 segments. Use "track.get_subtrack" to get subtracks.')
                end
            end
            
            obj.par = value;
        end
    end
    
    methods(Static)
        h_channel = get_los_channels( h_parset, precision, return_coeff, tx_array_mask )
    end
end

