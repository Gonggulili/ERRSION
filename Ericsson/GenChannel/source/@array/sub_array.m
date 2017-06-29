function a = sub_array( h_array, mask )
%SUB_ARRAY Generates a sub-array with the given array indices
%
%   This function creates a copy of the given array with only the selected
%   elements specified in mask.
%
% QuaDRiGa Copyright (C) 2011-2012 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
%
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

if ~isempty( setdiff( mask , 1:h_array.no_elements ) )
    error('The indices specified in mask do not exist in the array.')
end

h_array.uncompress;

a = array;

tmp = sprintf( '%d,',mask );
a.name =  [ h_array.name,'; El. ',tmp(1:end-1) ];

a.precision                 = h_array.precision;
a.interpolation_method      = h_array.interpolation_method;
a.elevation_grid            = h_array.elevation_grid;
a.azimuth_grid              = h_array.azimuth_grid;
a.no_elements               = numel( mask );
a.element_position          = h_array.element_position( :,mask );
a.Fa                        = h_array.Fa( :,:,mask );
a.Fb                        = h_array.Fb( :,:,mask );

if numel( h_array.Fc ) ~= 1
    a.Fc                    = h_array.Fc( :,:,mask );
end

if all( size( h_array.coupling ) == [1 1]*h_array.no_elements )
    a.coupling = h_array.coupling( mask,mask );
end

end
