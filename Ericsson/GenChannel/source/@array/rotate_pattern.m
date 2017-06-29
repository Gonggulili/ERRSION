function rotate_pattern( h_array, deg, rotaxis, element, usage )
%ROTATE_PATTERN Rotates antenna patterns
%
%   cp = ROTATE_PATTERN matches the pattern to the orientation vector. Use the
%   function "estimate_pol_vector" to determine best matching orientation
%   vectors first.
%
%   ROTATE_PATTERN( deg ) rotates the beam patterns of all elements around
%   the y-Axis in Cartesian coordinates. The value deg is the rotation
%   angle in degrees ranging from -180 to 180°.
%
%   ROTATE_PATTERN( deg,rotaxis ) also specifies the rotation axis x,y or z.
%   ROTATE_PATTERN( deg,rotaxis,element ) rotates only the given element.
%
%   ROTATE_PATTERN( deg,rotaxis,element,usage )
%   The optional parameter 'usage' can limit the rotation procedure either
%   to the pattern or polarization. Possible values are:
%
%       0: Rotate both (pattern+polarization) - default
%       1: Rotate only pattern
%       2: Rotate only polarization
%
%   Pattern rotation provides the option to assemble antenna arrays out of single
%   elements. By setting the element_position property of an array object, elements
%   can be placed at different coordinates. In order to freely design arbitrary
%   array configurations, however, elements often need to be rotated (e.g. to
%   assemble a +/- 45° cross-polarized array out of single dipoles). This
%   functionality is provided here.
%
% QuaDRiGa Copyright (C) 2011-2014 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
%
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% Parse input arguments
if exist('deg','var')
    if ~( all(size(deg) == [1 1]) && isnumeric(deg) ...
            && isreal(deg) )
        error('??? "deg" must be integer and real')
    end
else
    deg = 0;
end

if exist('rotaxis','var')
    if ~( ischar(rotaxis) && ...
            (strcmp(rotaxis,'x') || strcmp(rotaxis,'y')  || strcmp(rotaxis,'z')) )
        error('??? "rotaxis" can only be x,y or z.')
    end
else
    rotaxis = 'y';
end

if exist('element','var') && ~isempty( element )
    if ~( size(element,1) == 1 && isnumeric(element) ...
            &&  all( mod(element,1)==0 ) && min(element) > 0 && max(element)<=h_array.no_elements )
        error('??? "element" must be integer > 0 and can not exceed array size')
    end
else
    element = 1:h_array.no_elements;
end

if exist('usage','var')
    if ~( all(size(usage) == [1 1]) && isnumeric(usage) ...
            && any(usage == [0,1,2]) )
        error('??? "usage" must be 0,1 or 2')
    end
else
    usage = 0;
end

% Change the polarization basis to Polar-Spheric
polarization_basis = h_array.polarization_basis;
if ~strcmp( polarization_basis, 'polar-spheric' )
    h_array.change_pol_basis('polar-spheric');
end

% Get the angles.
phi   = h_array.azimuth_grid;
theta = h_array.elevation_grid';
no_az = h_array.no_az;
no_el = h_array.no_el;

% Rotation vectors are given in degree, but calculations are done in rad.
deg = deg/180*pi;

% Rotations are only allowed axis-wise where for each axis, another
% rotation matrix is given.
switch rotaxis
    case 'x'
        Rx = makehgtform('xrotate', deg);
        rot = Rx(1:3, 1:3);
    case 'y'
        Ry = makehgtform('yrotate', deg);
        rot = Ry(1:3, 1:3);
    case 'z'
        Rz = makehgtform('zrotate', deg);
        rot = Rz(1:3, 1:3);
end

% Calculate the transformation matrices for transforming the pattern from
% spherical coordinates to Cartesian coordinates.
B(1,:,:) = cos(theta)*cos(phi);             % Matrix-vector notation is faster
B(2,:,:) = cos(theta)*sin(phi);             % ... then meshgrid and sph2cart
B(3,:,:) = sin(theta)*ones(1,no_az);

A = rot.' * reshape(B, 3, []);
A = reshape( A.' , no_el,no_az,3 );

% Fix numeric bounds
A(A>1)  = 1;
A(A<-1) = -1;

% Calculate new angles
[phi_new, theta_new] = cart2sph( A(:,:,1), A(:,:,2), A(:,:,3) ) ;

% Angles might become complex, if the values in A are out of bound due
% to numeric offsets.
phi_new   = real( phi_new );
theta_new = real( theta_new );

% When the angles map to the poles of the pattern (i.e. theta = +/- 90
% degree), the pattern is not defined anymore. Here, we correct this by 
% using a small offset angle.
err_limit = 1e-5;

s = theta_new < -pi/2 + err_limit;
theta_new(s)  = -pi/2 + err_limit;

s = theta_new > pi/2 - err_limit;
theta_new(s)  = pi/2 - err_limit;

% Calculate the transformation for the polarization
if usage == 0 || usage == 2
    % Spherical basis vector in theta direction (original angles)
    Eth_o(1,:,:) = sin(theta)  * cos(phi);
    Eth_o(2,:,:) = sin(theta)  * sin(phi);
    Eth_o(3,:,:) = -cos(theta) * ones(1,no_az);
    Eth_o = reshape( Eth_o, 3, []);
    
    % Spherical basis vector in phi direction (original angles)
    Eph_o(1,:,:) = ones(no_el,1) * -sin(phi);
    Eph_o(2,:,:) = ones(no_el,1) * cos(phi);
    Eph_o(3,:,:) = zeros( no_el , no_az );
    Eph_o = reshape( Eph_o, 3, []);
    
    % Spherical basis vector in theta direction (new angles)
    Eth_n(1,:,:) = sin(theta_new) .* cos(phi_new);
    Eth_n(2,:,:) = sin(theta_new) .* sin(phi_new);
    Eth_n(3,:,:) = -cos(theta_new);
    Eth_n = reshape( Eth_n, 3, []);
    
    tmp = rot * Eth_n;
    cos_psi = sum( Eth_o .* tmp , 1 );
    sin_psi = sum( Eph_o .* tmp , 1 );
    
    cos_psi = reshape( cos_psi, no_el,no_az );
    sin_psi = reshape( sin_psi, no_el,no_az );
end

for n = 1:numel(element)
    
    % When we have the angles from the projection, we interpolate the
    % 3D field patterns to get the new (rotated) patterns. This is done
    % by the interpolation routine.
    
    % Interpolate the pattern
    if usage == 0 || usage == 1
        [ V,H ] = h_array.interpolate( phi_new, theta_new , element(n) );
    else
        V = h_array.Fa(:,:,element(n));
        H = h_array.Fb(:,:,element(n));
    end
    
    % Update the element position
    h_array.element_position(:,element(n)) = rot*h_array.element_position(:,element(n));
    
    % Transformation of the polarization
    if usage == 0 || usage == 2
        
        patV = cos_psi.*V - sin_psi.*H;
        patH = sin_psi.*V + cos_psi.*H;
        
        h_array.Fa(:,:,element(n)) = patV;
        h_array.Fb(:,:,element(n)) = patH;
        
    else % mode 1 - Rotate patterns only, but not the polarization
        
        h_array.Fa(:,:,element(n)) = V;
        h_array.Fb(:,:,element(n)) = H;
        
    end
end

% Change the polarization basis to Polar-Spheric
if ~strcmp( polarization_basis, 'polar-spheric' )
    h_array.change_pol_basis( polarization_basis );
end

end
