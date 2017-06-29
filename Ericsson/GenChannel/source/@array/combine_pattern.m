function combine_pattern( h_array , center_frequency )
%VIRTUAL_PATTERN Calculates a virtual pattern of the given array
%
%   When the inputs of an array are coupled (i.e. fed with the same
%   signal), then it is possible to combine the elements of the array. This
%   function calculates the virtual pattern in H/V by using the QuaDRiGa
%   simulator. 
%
%   Input:
%       center_frequency
%       The center frequency in [Hz]
%
% QuaDRiGa Copyright (C) 2011-2014 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
%
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

if ~exist('center_frequency','var')
    center_frequency = 2.6e9;           % Use default center frequency
    warning('QuaDRiGa:Array:virtual_pattern','Center Frequency is not given. Using 2.6 GHz.')
end

% The virtual pattern is calculated with respect to a reference antenna
% which we assume to be perfectly omni-directional with two polarizations.
% The reduction in resolution is made to save computing time.

ref = array('omni');
ref.set_grid( (-180:30:180)*pi/180 , (-90:30:90)*pi/180 );
ref.copy_element(1,2);
ref.Fb(:,:,2) = -ref.Fa(:,:,1);
ref.Fa(:,:,2) = 0;
ref.interpolation_method = 'nearest';
ref.precision = h_array.precision;

% Change the polarization basis to Polar-Spheric
polarization_basis = h_array.polarization_basis;
if ~strcmp( polarization_basis, 'polar-spheric' )
    h_array.change_pol_basis('polar-spheric');
end

% The input array is sampled at exactly the same angles for which the
% pattern is defined. Using nearest-neighbor interpolation thus speeds up
% the computations.
interpolation_method = h_array.interpolation_method;
h_array.interpolation_method = 'linear';

no_az = h_array.no_az;
no_el = h_array.no_el;

% The receiver positions are placed in 100 m distance in the same grid
% given by the elevation and azimuth angles in the original array.

phi   = h_array.azimuth_grid;
theta = h_array.elevation_grid';

B = zeros( 3,no_el,no_az );
B(1,:,:) = cos(theta)*cos(phi);
B(2,:,:) = cos(theta)*sin(phi);
B(3,:,:) = sin(theta)*ones(1,no_az);
B = 100*reshape(B, 3, []);

% Create a parameter_set object with the needed information. 
h_parset = parameter_set('LOSonly',[],0);
h_parset.simpar.center_frequency = center_frequency;
h_parset.tx_position = [0;0;0];
h_parset.no_positions = size(B,2);
h_parset.positions = B;
h_parset.rx_track = track('linear',0,0);
h_parset.rx_array = ref;
h_parset.tx_array = h_array;

% Obtain the channel coefficients
coeff = channel_builder.get_los_channels(h_parset,h_array.precision,'raw');
pat = permute( coeff , [3,2,1] ); % Map, Tx, Rx

% Write the output pattern
h_array.no_elements = size( pat,2 );
h_array.Fa = reshape( pat(:,:,1), no_el ,no_az , [] );
h_array.Fb = reshape( pat(:,:,2), no_el ,no_az , [] );
h_array.element_position = zeros(3,size( pat,2 ));
h_array.coupling = eye( h_array.no_elements );

% Restore old value
h_array.interpolation_method = interpolation_method;

% Change the polarization basis to Polar-Spheric
if ~strcmp( polarization_basis, 'polar-spheric' )
    h_array.change_pol_basis( polarization_basis );
end

end
