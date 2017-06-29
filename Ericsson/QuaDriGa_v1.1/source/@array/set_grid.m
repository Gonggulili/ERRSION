function set_grid( h_array , azimuth_grid , elevation_grid )
%SET_GRID Sets a new grid for azimuth and elevation and interpolates the pattern
%
%   set_grid( azimuth_grid , elevation_grid ) replaces the properties
%   "azimuth_grid" and "elevation_grid" of the antenna h_arrayect with the
%   values given here. Use this method when you also want to interpolate
%   the antenna patterns to the new grid. However, setting the class
%   parameters directly simply truncate the field patterns.
%
%   Example:
%   resolution = 3; 
%   a.set_grid( (-180:resolution:180)*pi/180 , (-90:resolution:90)*pi/180 );
%
% QuaDRiGa Copyright (C) 2011-2014 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
% 
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published 
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

if nargin < 3
   error('QuaDRiGa:Array:wrongNumberOfInputs',...
       'Wrong number of input arguments.');
end

if isempty( azimuth_grid ) || isempty( elevation_grid )
   error('QuaDRiGa:Array:wrongInputValue',...
       'Input arguments can not be empty.');
end

if ~( any( size(elevation_grid) == 1 ) && isnumeric(elevation_grid) && isreal(elevation_grid) &&...
        max(elevation_grid)<=pi/2 && min(elevation_grid)>=-pi/2 )
    error('QuaDRiGa:Array:wrongInputValue','??? "elevation_grid" must be a vector containing values between -pi/2 and pi/2')
end

if ~( any( size(azimuth_grid) == 1 ) && isnumeric(azimuth_grid) && isreal(azimuth_grid) &&...
        max(azimuth_grid)<=pi && min(azimuth_grid)>=-pi )
    error('QuaDRiGa:Array:wrongInputValue','??? "azimuth_grid" must be a vector containing values between -pi and pi')
end

el = repmat( elevation_grid' , 1 , numel(azimuth_grid) );
az = repmat( azimuth_grid , numel(elevation_grid) , 1 );

num = numel(h_array);
for n = 1 : num
    ind = h_array(n).eq( h_array ) ;
    ind( n ) = false;
    ind = find(ind,1);
    
    if ind < n
    else
        has_cartesian_basis = false;
        if strcmp( h_array(n).polarization_basis, 'cartesian' )
            has_cartesian_basis = false;
            h_array(n).change_pol_basis('polar-spheric');
        end
        if h_array(n).iscompressed
            h_array(n).uncompress;
        end
        
        iEl = elevation_grid <= max(h_array(n).elevation_grid) & ...
            elevation_grid >= min(h_array(n).elevation_grid);
        
        iAz = azimuth_grid <= max(h_array(n).azimuth_grid) & ...
            azimuth_grid >= min(h_array(n).azimuth_grid);
        
        [V,H] = h_array(n).interpolate( az(iEl,iAz) , el(iEl,iAz) , 1:h_array(n).no_elements  );
        h_array(n).elevation_grid = elevation_grid(iEl);
        h_array(n).azimuth_grid = azimuth_grid(iAz);
        h_array(n).Fa = V;
        h_array(n).Fb = H;
        
        if has_cartesian_basis
            h_array(n).change_pol_basis('cartesian');
        end
    end
end

end
