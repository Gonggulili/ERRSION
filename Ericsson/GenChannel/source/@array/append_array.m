function append_array( h_array, a )
%APPEND_ARRAY Appends an antenna array to the existing one
%
%   This method append the antenna array given in "a" to the existing
%   array object. If the polarization basis or the sampling grid do not
%   match, appropriate conversations or interpolations are performed.
%
% QuaDRiGa Copyright (C) 2011-2015 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
%
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

if ~isa( a, 'array' )
    error('QuaDRiGa:Array:append_array','"a" must be an array object.');
end

% Change the polarization basis to Polar-Spheric
polarization_basis = h_array.polarization_basis;
if ~strcmp( polarization_basis, 'polar-spheric' )
    h_array.change_pol_basis('polar-spheric');
end

if ~strcmp( a.polarization_basis, 'polar-spheric' )
    b = a.copy;
    b.change_pol_basis('polar-spheric');
else
    b = a;
end

% Match the sampling grids
if numel( b.elevation_grid ) ~= numel( h_array.elevation_grid ) || ...
        numel( b.azimuth_grid ) ~= numel( h_array.azimuth_grid ) || ...
        any( h_array.elevation_grid - b.elevation_grid > 1e-12 ) || ...
        any( h_array.azimuth_grid - b.azimuth_grid > 1e-12 )
    
    warning('QuaDRiGa:Array:append_array','The array sampling grids do not match, therefore the sampling grid of ''a'' will be changed and the pattern will be interpolated accordingly before appending.');
    b = b.copy;
    b.set_grid( h_array.azimuth_grid , h_array.elevation_grid );
end

cpl_in  = h_array.coupling;
cpl_app = a.coupling;

% Set the number of elements in the new array
i_el_in = h_array.no_elements;
h_array.no_elements = i_el_in + b.no_elements;

% Copy the pattern information
for i_el = 1 : a.no_elements
    i_el_out = i_el_in+i_el;
    h_array.Fa( : , : , i_el_out ) = b.Fa( :,:,i_el );
    h_array.Fb( : , : , i_el_out ) = b.Fb( :,:,i_el );
    h_array.element_position( :,i_el_out ) = b.element_position( :,i_el );
end

% Update the coupling matrix
nc = size(cpl_in);
ne = size(cpl_app);
cpl = zeros( nc(1)+ne(1) , nc(2)+ne(2));
cpl( 1:nc(1) , 1:nc(2) ) = cpl_in;
cpl( nc(1)+1:end , nc(2)+1:end ) = cpl_app;

h_array.coupling = cpl;

% Change the polarization basis to Polar-Spheric
if ~strcmp( polarization_basis, 'polar-spheric' )
    h_array.change_pol_basis( polarization_basis );
end

end
