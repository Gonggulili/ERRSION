function coupling = generate_multi( h_array, no_elements, spacing, tilt, optimize )
%ELECTRIC_TILT Generates a multi-element array with electric tilt
%
%   This function generates a vertically stacked multi-element array of the
%   given antenna object. The element spacing is relative to the wavelength
%   and an additional electric tilt can be applied. The output is stored in
%   the current antenna object.
%
%   For this method to work, you need to define a single antenna element
%   using "array.generate". Then, you can call "array.generate_multi" to
%   transform this element into a stacked multi-element array. The provided
%   antenna object can only have one element. The method returns an error
%   if a multi-element array is given.
%
%   Input:
%       no_elements
%       The number of "virtual" antenna elements stacked in elevation (z)
%       direction.
%
%       spacing
%       The element-spacing as a factor of the wavelength.
%       Default value: 0.5
%
%       tilt
%       An additional electric downtilt value in [deg].
%       Default value: 0
%
%       optimize
%       If this parameter is set to 1, the optimal beamformer is calculated.
%       Otherwise, the phases are calculated using geometric settings.
%       Default: 0 (geometric)
%
% QuaDRiGa Copyright (C) 2011-2015 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
%
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.


if ~exist('no_elements','var')
    error('QuaDRiGa:Array:generate_multi','Number of elements is not specified.');
elseif ~isreal( no_elements ) || no_elements < 2
    error('QuaDRiGa:Array:generate_multi','Number of elements must be larger than 1.');
end

s = simulation_parameters;
if ~exist('spacing','var') || isempty( spacing )
    spacing = 0.5 * s.wavelength;
    warning('QuaDRiGa:Array:generate_multi','Spacing is not given. Using 0.5*lambda.');
elseif ~isreal( spacing ) || spacing <= 0
    error('QuaDRiGa:Array:generate_multi','Spacing must be larger than 0.');
else
    spacing = spacing * s.wavelength;
end

if ~exist('tilt','var') || isempty( tilt )
    tilt = 0;
    warning('QuaDRiGa:Array:generate_multi','Tilt is not given. Using 0 degree.');
else
    % Convert tilt to [rad]
    tilt = tilt * pi/180;  
end

if ~exist('optimize','var') || isempty( optimize )
    optimize = 0;
end

% Change the polarization basis to Polar-Spheric
polarization_basis = h_array.polarization_basis;
if ~strcmp( polarization_basis, 'polar-spheric' )
    h_array.change_pol_basis('polar-spheric');
end

% Copy the basic elements
no_el_in = h_array.no_elements;
for n = 2:no_el_in
    h_array.copy_element( n, (n-1)*no_elements+1 );
end

% Set the element spacing
el_pos = (0:no_elements-1) * spacing;
el_pos = el_pos - mean(el_pos);
for n = 1:no_el_in
    ind = (n-1)*no_elements + (1:no_elements);
    h_array.copy_element( (n-1)*no_elements+1 , ind );
    h_array.element_position(3,ind) = el_pos;
end


if optimize
    % Use the channel model to calculate the optimal weights.
    
    % calculate the phases
    h_parset = parameter_set('LOSonly');
    
    % Create h_array dual-polarized Rx-Array
    h_parset.tx_array = h_array;
    h_parset.rx_array.Fa = 1/sqrt(2) * h_parset.rx_array.Fa;
    h_parset.rx_array.Fb = h_parset.rx_array.Fa;
    
    h_parset.tx_position = [0,0,0]';
    
    % The virtual receiver position
    tmp = exp(1j*tilt) * 1000;
    h_parset.positions(1) = real( tmp );
    h_parset.positions(3) = imag( tmp );
    h_parset.rx_track.generate('linear',0,0);
    
    % Get the phases
    h_channel = channel_builder.get_los_channels(h_parset);
    
    coupling  = exp( 1j*angle( h_channel.coeff )' );
    coupling  = reshape( coupling, no_elements , no_el_in );
    coupling  = coupling * sqrt(1/no_elements);
    
else
    % Geometric calculation
    
    C = -2*pi*sin(tilt) * h_array.element_position(3,:) / s.wavelength;
    C = exp( 1j*C );
    C = reshape( C, no_elements , no_el_in );
    coupling = C * sqrt(1/no_elements);
    
end

coupling_array = zeros( no_el_in*no_elements, no_el_in );
for n = 1:no_el_in
    ind = (n-1)*no_elements + (1:no_elements);
    coupling_array(ind,n) = coupling(:,n);
end

h_array.coupling = coupling_array;

% calculate the effective pattern
h_array.combine_pattern( s.center_frequency );

% Change the polarization basis to Polar-Spheric
if ~strcmp( polarization_basis, 'polar-spheric' )
    h_array.change_pol_basis( polarization_basis );
end

end
