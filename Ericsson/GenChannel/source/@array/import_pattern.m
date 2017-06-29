function h_array = import_pattern( fVi, fHi , azimuth_grid , elevation_grid )
%IMPORT_PATTERN Converts antenna field patterns into a QuaDRiGa array object
%
%   This function converts any antenna field pattern into a QuaDRiGa
%   antenna array object.
%
%   The input variables are:
%
%    fVi
%      The field pattern(s) for the vertical polarization given in
%      spherical coordinates. The first dimension corresponds to the
%      elevation angle (ranging from -90 to 90 degrees). The second
%      dimension is for the azimuth angle (ranging from -180 to 180
%      degrees). The third dimension belongs to the element number. The
%      default resolution is 1 degree - hence, the default size of fVi is
%      <181x361x1>. If a different resolution is given, the optional
%      variables "azimuth_grid" and "elevation_grid" must be defined.
%
%    fHi
%      The field pattern(s) for the horizontal polarization given in
%      spherical coordinates. "fHi" can be empty if no horizontal response
%      is given. If it is given, then "fHi" must have the same size as
%      "fVi".
%
%    azimuth_grid
%      A vector specifying the azimuth sampling points of the patterns in
%      units of radians (raging from -pi to pi). This value only needs to
%      be defined if the patterns do not have the default size.
%
%    elevation_grid
%      A vector specifying the elevation sampling points of the patterns in
%      units of radians (raging from -pi/2 to pi/2). This value only needs
%      to be defined if the patterns do not have the default size.
%
%
%   The output variables are:
%
%     h_array
%       The QuaDRiGa antenna array object generated from the field patterns.
%
%
% QuaDRiGa Copyright (C) 2011-2014 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
%
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% Parse input: fVi
if ~exist('fVi','var') || isempty(fVi)
    error('Vertical field pattern "fVi" is not defined.')
end

% Parse input: fHi
if exist('fHi','var') && ~isempty(fHi)
    if any( size(fVi) ~= size(fHi) )
        error('Horizontal field pattern "fHi" has different size than vertical pattern "fVi".')
    end
else
    fHi = zeros( size(fVi) );
end

% Parse input: azimuth_grid
if exist('azimuth_grid','var') && ~isempty(azimuth_grid)
    if ~( any( size(azimuth_grid) == 1 ) && isnumeric(azimuth_grid) && isreal(azimuth_grid) &&...
            max(azimuth_grid)<=pi && min(azimuth_grid)>=-pi )
        error('??? "azimuth_grid" must be a vector containing values between -pi and pi')
    end
elseif size(fVi,2) ~= 361
    error('??? "azimuth_grid" undefined');
else
    azimuth_grid = (-180:180)*pi/180;
end

% Parse input: elevation_grid
if exist('elevation_grid','var') && ~isempty(elevation_grid)
    if ~( any( size(elevation_grid) == 1 ) && isnumeric(elevation_grid) && isreal(elevation_grid) &&...
            max(elevation_grid)<=pi/2 && min(elevation_grid)>=-pi/2 )
        error('??? "elevation_grid" must be a vector containing values between -pi/2 and pi/2')
    end
elseif size(fVi,1) ~= 181
    error('??? "elevation_grid" undefined');
else
    elevation_grid = (-90:90)*pi/180;
end

no_input_elements = size(fVi,3);

h_array = array;
h_array.no_elements = no_input_elements;
h_array.azimuth_grid = azimuth_grid;
h_array.elevation_grid = elevation_grid;
h_array.Fa = fVi;
h_array.Fb = fHi;

end
