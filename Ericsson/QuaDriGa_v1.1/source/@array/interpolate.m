function [ V, H, dist] = interpolate( h_array, azimuth, elevation, element,...
    azimuth_grid, elevation_grid, Fa, Fb, element_position )
%INTERPOLATE Interpolates the field pattern
%
%   [V,H] = INTERPOLATE( azimuth,elevation ) Interpolates the field pattern to the
%   given azimuth and elevation values. Azimuth and elevation are angles in the
%   spheric coordinate system. Azimuth can range from -Pi to Pi and elevation can
%   have values in between -Pi/2 and Pi/2. Values outside this range will be
%   converted accordingly. Important: All angles must be given in radians!
%
%   [V,H] = INTERPOLATE( azimuth,elevation,element ) also specifies the element
%   range for which the pattern will be interpolated.
%
%   [V,H,CP] = INTERPOLATE( azimuth,elevation,element ) also returns the
%   interpolated phase map of the antenna.
%
%   [V,H,CP,dist] = INTERPOLATE( azimuth,elevation,element ) also returns the
%   effective distances between the antenna elements when seen from the
%   direction of the incident path. The distance is calculated by an
%   projection of the array positions on the normal plane of the incident
%   path.
%
%   Note: Interpolation of the beam patterns is very (!!!) computing intensive. It
%   must be performed several thousands of times during a simulation run. It is
%   highly recommended to use linear interpolation for this task since
%   this method is optimized for speed. Spline interpolation calls the
%   MATLAB internal interpolation function which is more than 10 times slower.
%   To enable linear interpolation, set the 'interpolation_method'
%   property of the array object to 'linear'.
%
% QuaDRiGa Copyright (C) 2011-2015 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
%
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% Parse input arguments
if nargin < 5
    if nargin < 3
        error('??? You need to specify azimuth and elevation for the pattern interpolation.');
    elseif nargin < 4
        element = 1:h_array.no_elements;
    end
    
    if strcmp( h_array.polarization_basis, 'cartesian' )
        error('??? Array interpolation does not work for Cartesian polarization basis.');
    end
    
    azimuth_grid  = h_array.azimuth_grid;
    elevation_grid = h_array.elevation_grid;
    Fa = h_array.Fa;
    Fb = h_array.Fb;
    element_position = h_array.element_position(:,element);
end

if isa( Fa , 'single' )
    single_precision = true;
else
    single_precision = false;
end

% Get initial values
dims = size( azimuth );
no_values = numel( azimuth );
no_element = numel( element );
no_az = h_array.no_az;
no_el = h_array.no_el;

% Reduce dimensions of arrays
azimuth = reshape( azimuth,1,no_values );
elevation = reshape( elevation,1,no_values );

% elevation may be outside of (-pi/2,pi/2).
% Mapping of angles to (-pi/2,pi/2) and changing the azimuthal orientation of rays
% with original elevation angles within +-(pi/2,3pi/2)+n*2*pi. 29.9.2008 EsaK

ind = floor( elevation/pi + 0.5 );
elevation = (-1).^(ind + 1).*(pi*ind - elevation);

ind = logical( mod(ind,2) );
azimuth(ind) = azimuth(ind) + pi;

azimuth = mod( azimuth+pi , 2*pi )-pi;

% Initialize variables
if single_precision
    V = zeros( no_values ,no_element, 'single' );
    H = zeros( no_values ,no_element, 'single' );
else
    V = zeros( no_values ,no_element );
    H = zeros( no_values ,no_element );
end

% Interpolate field patterns
switch h_array.interpolation_method
    case 'nearest'
        
        % Determine the nearest location of xi in x
        [tmp,b] = sort( azimuth );
        [~,a]   = sort( [azimuth_grid,tmp] );
        ui      = 1:(no_az + no_values);
        ui(a)   = ui;
        ui      = ui(no_az+1:end) - (1:no_values);
        ui(b)   = ui;
        
        % Determine the nearest location of yi in y
        [tmp,b] = sort( elevation );
        [~,a]   = sort( [elevation_grid,tmp] );
        vi      = 1:(no_el + no_values);
        vi(a)   = vi;
        vi      = vi(no_el+1:end) - (1:no_values);
        vi(b)   = vi;
        
        % Determine the index of the element and get the corresponding pattern entry
        vi = vi+(ui-1)*no_el;
        for n=1:no_element
            ndx = (element(n)-1)*no_az*no_el + vi;
            V(:,n) = Fa(ndx);
            H(:,n) = Fb(ndx);
        end
        
    case 'linear'
        tmp_1_no_vals = uint32( (1:no_values) );
        
        % Determine the nearest location of xi in x and the difference to
        % the next point
        [tmp,b] = sort( azimuth );
        [~,a]   = sort( [azimuth_grid,tmp] );
        ui      = uint32( 1:(no_az + no_values) );
        ui(a)   = ui;
        ui      = ui(no_az+1:end) - tmp_1_no_vals;
        ui(b)   = ui;
        ui( ui==no_az ) = no_az-1;
        ui( ui==0 ) = 1;
        uin     = ui+1;
        u       = (azimuth-azimuth_grid(ui))./( azimuth_grid(uin)-azimuth_grid(ui) );
        u       = u';
        
        % Determine the nearest location of yi in y and the difference to
        % the next point
        [tmp,b] = sort( elevation );
        [~,a]   = sort( [elevation_grid,tmp] );
        vi      = uint32( 1:(no_el + no_values) );
        vi(a)   = vi;
        vi      = vi(no_el+1:end) - tmp_1_no_vals;
        vi(b)   = vi;
        vi( vi==no_el ) = no_el-1;
        vi( vi==0 ) = 1;
        vin     = vi+1;
        v       = (elevation-elevation_grid(vi))./( elevation_grid(vin)-elevation_grid(vi) );
        v       = v';
        
        % Calculate the scaling coefficients
        c1 = (1-v).*(1-u);
        c2 = (1-v).*u;
        c3 = v.*(1-u);
        c4 = v.*u;
        
        % Determine the indices of the elements
        pa = vi  + ( ui  -1 )*no_el;
        pb = vi  + ( uin -1 )*no_el;
        pc = vin + ( ui  -1 )*no_el;
        pd = vin + ( uin -1 )*no_el;
        
        pX = [pa,pb,pc,pd].';
        pY = uint32( (element-1)*no_az*no_el );
        
        tr = true( no_values,1 );
        fl = false( no_values,1 );
        i1 = [tr;fl;fl;fl];
        i2 = [fl;tr;fl;fl];
        i3 = [fl;fl;tr;fl];
        i4 = [fl;fl;fl;tr];
        
        % Interpolate
        for n = 1 : no_element
            ndx = pY(n) + pX;
            
            a = Fa( ndx );
            V(:,n) = c1.*a(i1) + c2.*a(i2) + c3.*a(i3) + c4.*a(i4);
            
            a = Fb( ndx );
            H(:,n) = c1.*a(i1) + c2.*a(i2) + c3.*a(i3) + c4.*a(i4);
        end
        
    otherwise
        
        % Expand the interpolation grid to the left and right to reduce errors at
        % pi.
        
        a = numel(azimuth_grid);
        ind = [ a-5:a-1 1:a 2:5  ];
        agn = [ azimuth_grid(a-5:a-1)-2*pi azimuth_grid azimuth_grid(2:5)+2*pi ];
        
        % Call the MATLAB-internal interpolation functions
        for n = 1:no_element
            V(:,n) = interp2( agn , elevation_grid ,...
                Fa(:,ind,element(n)) , azimuth , elevation , h_array.interpolation_method );
            
            H(:,n) = interp2( agn , elevation_grid ,...
                Fb(:,ind,element(n)) , azimuth , elevation , h_array.interpolation_method );
        end
end

% Remap output to match input dimensions

V = reshape( V, [dims,no_element] );
H = reshape( H, [dims,no_element] );

% The effective distances are only needed when we calculate the channel
% coefficients. Therefore, it is also provided by the interpolation
% function.

if nargout > 2
    % Compute the effective distances for each path
    
    % Transform incoming angles from spherical coordinates to Cartesian
    % coordinates.
    if any(any( element_position~=0 ))
        if single_precision
            A = zeros(3,no_values,'single');
        else
            A = zeros(3,no_values);
        end
        [A(1, :), A(2, :), A(3, :)] = sph2cart( azimuth, elevation, 1 );
    end
    
    % The distances are calculated by a parallel projection.
    % See: http://de.wikipedia.org/wiki/Parallelprojektion
    % The normal vector of the projection plane is given in X.
    % The original point P is the element position defined by the array
    % geometry. The origin of the projection plane is the origin of the
    % array coordinate system. The projection result is B.
    
    if single_precision
        dist = zeros( no_values, no_element,'single' );
    else
        dist = zeros( no_values, no_element );
    end
    
    for n = 1 : no_element
        E = element_position(:,n);
        
        if any( E ~= 0 )            % The distance is 0 if P is in the origin
            tmp = E'*A;             % cos(E,A)
            D = tmp([1 1 1],:) .* A;
            dist(:,n) = -sign( A(1,:).*D(1,:) + A(2,:).*D(2,:) + A(3,:).*D(3,:) ) .* ...
                sqrt(sum( D.^2 ));
            
            % This is the same as
            % for m = 1:no_values
            %    D = -( -E.' * A(:,m) ) * A(:,m) ;
            %    dist(m,n) = sign(D'*A(:,m)) * sqrt(sum(  D.^2 ));
            % end
        end
    end
    
    % Return the output
    dist = reshape( dist, [dims,no_element] );
end

end
