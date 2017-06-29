function par = generate(h_array, array_type, element, Ain, Bin, Cin, Din, Ein, Fin )
%GENERATE Generates predefined arrays
%
%   GENERATE( array_type ) generates a new, predefined antenna array given by
%   array_type. Currently supported are:
%
%   'omni'				An isotropic radiator with vertical polarization.
%
%   'dipole'/'short-dipole'		A short dipole radiating with vertical polarization.
%
%   'half-wave-dipole'	A half-wave dipole radiating with vertical polarization.
%
%   'patch'  			An patch antenna with 90 deg opening in azimuth and
%           			elevation.
%
%   'custom'			An antenna with a custom gain in elevation and azimuth.
%                       Ain = 3dB beam width in azimuth direction
%                       Bin = 3dB beam width in elevation direction
%                       Cin = Isotropic gain (linear scale) at the back of the antenna
%                       The values A,B,C and D for the parametric antenna are returned.
%
%   '3gpp-macro'        An antenna with a custom gain in elevation and azimuth.
%                       Ain = Half-Power in azimuth direction (default = 70 deg)
%                       Bin = Half-Power in elevation direction (default = 10 deg)
%                       Cin = Front-to-back ratio (default = 25 dB)
%                       Din = Electrical downtilt (default = 15 deg)
%                       See 3GPP TR 36.814 V9.0.0 (2010-03), Table
%                       A.2.1.1-2, Page 59
%
%   '3gpp-3d'           The antenna model for the 3GPP-3D channel model.
%                       See 3GPP TR 36.873, v12.2.0, page 16
%                       Ain = Number of vertical elements (M)
%                       Bin = Number of horizontal elements (N)
%                       Cin = The center frequency in [Hz]
%                       Din = Polarization indicator
%                           1: K=1, vertical polarization only
%                           2: K=1, H/V polarized elements
%                           3: K=1, +/-45 degree polarized elements
%                           4: K=M, vertical polarization only
%                           5: K=M, H/V polarized elements
%                           6: K=M, +/-45 degree polarized elements
%                       Ein = The electric tile angle in [deg] for Cin = {4,5,6}
%                       Fin = Element spacing in [lambda], Default: 0.5
%
%   'parametric'        An antenna following the function
%                       F = A * sqrt( B + (1-B)*( cos(theta) )^C * exp( -D
%                       * phi^2 ) );
%
%   'xpol'   			Two elements with ideal isotropic patterns (vertical
%           			polarization). The second element is tilted by 90 degree.
%
%   'rhcp-dipole'   	Two crossed dipoles with one port. The signal on the
%           			second element (horizontal) is shifted by -90 deg out of phase.
%           			The two elements thus create a RHCP signal.
%
%   'lhcp-dipole'   	Two crossed dipoles with one port. The signal on the
%           			second element (horizontal) is shifted by +90 deg out of phase.
%           			The two elements thus create a LHCP signal.
%
%   'lhcp-rhcp-dipole'  Two crossed dipoles. For input port 1, the signal on
%           			the second element (horizontal) is shifted by +90 deg out of
%           			phase. For input port 2, the the signal on the second element
%           			(horizontal) is shifted by -90 deg out of phase. Port 1 thus
%           			transmits a LHCP signal and port 2 transmits a RHCP signal.
%
%   'ulaX'    			Unified linear arrays composed of omni-antennas (vertical
%           			polarization)  with 10 cm element distance. X can be 2,4 or 8.
%
% QuaDRiGa Copyright (C) 2011-2016 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
%
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

supported_types = array.supported_types;
if ~exist( 'array_type' , 'var' ) || isempty(array_type)
    array_type = 'omni';
elseif ~( ischar(array_type) && any( strcmpi(array_type,supported_types)) )
    str = ['Array type "',array_type,'" not found. Supported types are: '];
    no = numel(supported_types);
    for n = 1:no
        str = [str,supported_types{n}];
        if n<no
            str = [str,', '];
        end
    end
    error(str);
end

if ~exist( 'element' , 'var' ) || isempty(element)
    element = 0;
else
    if ~( size(element,1) == 1 && isnumeric(element) && all(size(element) == [1 1]) ...
            &&  all( mod(element,1)==0 ) && min(element) >= 0 )
        error('??? "element" must be scalar, integer >= 0 and can not exceed array size')
    end
end

par = [];

if numel(h_array) > 1
    % Do for parameter_set_arrays
    for n=1:numel(h_array)
        h_array(n).generate( array_type,1 );
    end
    
else
    if element > h_array.no_elements
        h_array.no_elements = element;
    end
    
    switch array_type
        
        case 'omni'
            if element == 0 || isempty(h_array.no_elements)
                h_array.no_elements                 = 1;
                h_array.elevation_grid              = (-90:90)*pi/180;
                h_array.azimuth_grid                = (-180:180)*pi/180;
                h_array.element_position            = zeros( 3,1 );
                h_array.Fa      = ones( 181,361 );
                h_array.Fb    = zeros( 181,361 );
                h_array.coupling                    = 1;
            else
                h_array.element_position(:,element)           = zeros( 3,1 );
                h_array.Fa(:,:,element)   = ones( h_array.no_el,h_array.no_az );
                h_array.Fb(:,:,element) = zeros( h_array.no_el,h_array.no_az );
            end
            
        case {'short-dipole', 'dipole'}
            if element == 0 || isempty(h_array.no_elements)
                h_array.generate('omni');
                h_array.interpolation_method        = 'linear';
                e = 1;
            else
                e = element;
            end
            [~, theta_grid] = meshgrid(h_array.azimuth_grid, h_array.elevation_grid);
            
            % Short dipole
            E_theta = cos( (1 - 1e-6) * theta_grid );
            E_phi = zeros(size(E_theta));
            
            P = E_theta.^2 + E_phi.^2;      % Calculate radiation power pattern
            P_max = max(max(P));            % Normalize by max value
            P = P ./ P_max;
            
            % Calculate the Gain
            gain_lin = sum(sum( cos(theta_grid) )) / sum(sum( P.*cos(theta_grid) ));
            
            % Normalize by Gain
            E_theta = E_theta .* sqrt(gain_lin./P_max);
            E_phi = E_phi .* sqrt(gain_lin./P_max);
            
            h_array.Fa(:,:,e) = E_theta;
            h_array.Fb(:,:,e) = E_phi;
            
        case 'half-wave-dipole'
            if element == 0 || isempty(h_array.no_elements)
                h_array.generate('omni');
                h_array.interpolation_method        = 'linear';
                e = 1;
            else
                e = element;
            end
            [~, theta_grid] = meshgrid(h_array.azimuth_grid, h_array.elevation_grid);
            
            % Half-Wave dipole
            E_theta = cos( pi/2*sin((1 - 1e-6) * theta_grid )) ./ cos((1 - 1e-6) * theta_grid);
            E_theta( isnan(E_theta) ) = 0;
            E_phi = zeros(size(E_theta));
            
            P = E_theta.^2 + E_phi.^2;      % Calculate radiation power pattern
            P_max = max(max(P));            % Normalize by max value
            P = P ./ P_max;
            
            % Calculate the Gain
            gain_lin = sum(sum( cos(theta_grid) )) / sum(sum( P.*cos(theta_grid) ));
            
            % Normalize by Gain
            E_theta = E_theta .* sqrt(gain_lin./P_max);
            E_phi = E_phi .* sqrt(gain_lin./P_max);
            
            h_array.Fa(:,:,e) = E_theta;
            h_array.Fb(:,:,e) = E_phi;
            
        case 'patch'
            if element == 0 || isempty(h_array.no_elements)
                h_array.generate('custom',0,90,90,0);
            else
                h_array.generate('custom',element,90,90,0);
            end
            
        case 'custom'
            
            % Set input variables
            if nargin < 4
                phi_3dB = 120;
            else
                phi_3dB = Ain;
            end
            
            if nargin < 5
                theta_3dB = 120;
            else
                theta_3dB = Bin;
            end
            
            if nargin < 6
                rear_gain = 0.1;
            else
                rear_gain = Cin;
            end
            
            if ~( size(phi_3dB,1) == 1 && isnumeric(phi_3dB) && isreal(phi_3dB) &&...
                    all(size(phi_3dB) == [1 1]) )
                error('Azimuth HPBW (phi_3dB) has invalid value.')
            end
            
            if ~( size(theta_3dB,1) == 1 && isnumeric(theta_3dB) && isreal(theta_3dB) &&...
                    all(size(theta_3dB) == [1 1]) )
                error('Elevation HPBW (theta_3dB) has invalid value.')
            end
            
            if ~( size(rear_gain,1) == 1 && isnumeric(rear_gain) && isreal(rear_gain) && ...
                    all(size(rear_gain) == [1 1]) && rear_gain>=0 && rear_gain<0.5)
                error('Front-to-back ratio (rear_gain) has invalid value.')
            end
            
            if element == 0 || isempty(h_array.no_elements)
                h_array.generate('omni');
                h_array.interpolation_method        = 'linear';
                e = 1;
            else
                e = element;
            end
            
            par.A = 0;
            par.B = 0;
            par.C = 0;
            par.D = 0;
            
            % Calculate the azimuth response
            phi = h_array.azimuth_grid;
            ind = find(phi/pi*180 >= phi_3dB/2, 1);
            
            a   = 1;        % Initial angle
            dm  = 0.5;      % Step size
            x   = inf;
            delta = Inf;
            ddir = +1;
            lp = 1;
            while lp < 5000 && delta > 1e-7
                if lp > 1
                    an = a + ddir * dm;
                    delta = abs(a-an);
                else
                    an = a;
                end
                
                C = rear_gain + (1 - rear_gain) * exp(-an * phi.^2);
                xn = abs(C(ind) - 0.5);
                
                if xn < x
                    a = an;
                    x = xn;
                else
                    ddir = -ddir;
                    dm = 0.382 * dm;
                end
                lp = lp + 1;
            end
            C = exp(-an * phi.^2);
            par.D = an;
            
            % Calculate the elevation response
            theta = h_array.elevation_grid;
            ind = find(theta/pi*180 >= theta_3dB/2, 1);
            
            a   = 1;        % Initial angle
            dm  = 0.5;      % Step size
            x   = inf;
            delta = Inf;
            ddir = +1;
            lp = 1;
            while lp < 5000 && delta > 1e-7
                if lp > 1;
                    an = a + ddir * dm;
                    delta = abs(a-an);
                else
                    an = a;
                end
                
                D = cos(theta).^an;
                xn = abs(D(ind) - 0.5);
                
                if xn < x
                    a = an;
                    x = xn;
                else
                    ddir = -ddir;
                    dm = 0.382 * dm;
                end
                lp = lp + 1;
            end
            D = cos(theta).^an;
            par.C = an;
            
            par.B = rear_gain;
            
            P = zeros(181,361);
            for a = 1:181
                for b = 1:361
                    P(a, b) = D(a) * C(b);
                end
            end
            P = rear_gain + (1-rear_gain)*P;
            
            E_theta =  sqrt(P);
            
            [~, theta_grid] = meshgrid(h_array.azimuth_grid, h_array.elevation_grid);
            
            P = E_theta.^2;         % Calculate radiation power pattern
            P_max = max(max(P));    % Normalize by max value
            P = P ./ P_max;
            
            % Calculate the Gain
            gain_lin = sum(sum( cos(theta_grid) )) / sum(sum( P.*cos(theta_grid) ));
            par.A = sqrt(gain_lin./P_max);
            
            E_theta = E_theta .* sqrt(gain_lin./P_max);
            
            h_array.Fa(:,:,e) = E_theta;
            h_array.Fb(:,:,e) = zeros(h_array.no_el, h_array.no_az);
            
        case '3gpp-macro'
            % Set input variables
            if nargin < 4
                phi_3dB = 70;
            else
                phi_3dB = Ain;
            end
            
            if nargin < 5
                theta_3dB = 10;
            else
                theta_3dB = Bin;
            end
            
            if nargin < 6
                rear_gain = 25;
            else
                rear_gain = Cin;
            end
            
            if nargin < 7
                electric_tilt = 15;
            else
                electric_tilt = Din;
            end
            
            if ~( size(phi_3dB,1) == 1 && isnumeric(phi_3dB) && isreal(phi_3dB) &&...
                    all(size(phi_3dB) == [1 1]) )
                error('Azimuth HPBW (phi_3dB) has invalid value.')
            end
            
            if ~( size(theta_3dB,1) == 1 && isnumeric(theta_3dB) && isreal(theta_3dB) &&...
                    all(size(theta_3dB) == [1 1]) )
                error('Elevation HPBW (theta_3dB) has invalid value.')
            end
            
            if ~( size(rear_gain,1) == 1 && isnumeric(theta_3dB) && isreal(rear_gain) && ...
                    all(size(rear_gain) == [1 1]) && rear_gain>=0 )
                error('Front-to-back ratio (rear_gain) has invalid value.')
            end
            
            if ~( size(electric_tilt,1) == 1 && isnumeric(electric_tilt) && isreal(electric_tilt) &&...
                    all(size(electric_tilt) == [1 1]) )
                error('Electric tilt has invalid value.')
            end
            
            if element == 0 || isempty(h_array.no_elements)
                h_array.generate('omni');
                h_array.interpolation_method = 'linear';
                e = 1;
            else
                e = element;
            end
            
            phi = h_array.azimuth_grid*180/pi;
            Ah  = -min( 12*(phi./phi_3dB).^2 , rear_gain );
            
            theta = h_array.elevation_grid.'*180/pi;
            Av  = -min( 12*((theta+electric_tilt)./theta_3dB).^2 , rear_gain-5 );
            
            A = -min( -Ah(ones(h_array.no_el,1),:) - Av(:,ones(h_array.no_az,1)) , rear_gain );
            
            h_array.Fa(:,:,e) = sqrt( 10.^(0.1*A) );
            h_array.normalize_gain(e);
            
        case '3gpp-3d'
            
            % Set inputs
            if ~exist('Ain','var') || isempty( Ain )
                M = 10;
            else
                M = Ain;
            end
            if ~exist('Bin','var') || isempty( Bin )
                N = 10;
            else
                N = Bin;
            end
            if ~exist('Cin','var') || isempty( Cin )
                center_freq = 2e9;
            else
                center_freq = Cin;
            end
            if ~exist('Din','var') || isempty( Din )
                pol = 1;
            else
                pol = Din;
            end
            if ~exist('Ein','var') || isempty( Ein )
                tilt = 0;
            else
                tilt = Ein;
            end
            if ~exist('Fin','var') || isempty( Fin )
                spacing = 0.5;
            else
                spacing = Fin;
            end
            
            % Create single element pattern
            h_array.generate('omni');
            
            % Antenna element vertical radiation pattern (dB)
            theta = h_array.elevation_grid.'*180/pi;
            Av  = -min( 12*(theta./65).^2 , 30 );
            
            % Antenna element horizontal radiation pattern (dB)
            phi = h_array.azimuth_grid*180/pi;
            Ah  = -min( 12*(phi./65).^2 , 30 );
            
            % Combining method for 3D antenna element pattern (dB)
            A = -min( -Ah(ones(h_array.no_el,1),:) - Av(:,ones(h_array.no_az,1)) , 30 );
            
            % Set pattern
            h_array.Fa = sqrt( 10.^(0.1*A) );
            
            % Maximum directional gain of an antenna element is 8 dB
            h_array.normalize_gain(1,8);
            
            % Polarization
            switch pol
                case {2,5} % H / V polarization (0, 90 deg slant)
                    h_array.copy_element(1,2);
                    h_array.rotate_pattern(90,'x',2,2);
                    
                case {3,6} % +/- 45 deg polarization
                    h_array.copy_element(1,2);
                    h_array.rotate_pattern(45,'x',1,2);
                    h_array.rotate_pattern(-45,'x',2,2);
            end
            
            % Coupling of vertically stacked elements
            if pol >= 4
                w = h_array.generate_multi( M, spacing, tilt );
                M = 1;
            end
            
            if N > 1 || M > 1
          % Calculate the wavelength
                s = simulation_parameters;
                s.center_frequency = center_freq;
                lambda = s.wavelength;
                
                % Copy elements
                T = h_array.no_elements;
                h_array.no_elements = T*N*M;
                for t=2:T
                    ii = T+t : T : T*N*M;
                    ij = ones(1,numel(ii))*t;
                    h_array.Fa(:,:,ii) = h_array.Fa(:,:,ij);
                    h_array.Fb(:,:,ii) = h_array.Fb(:,:,ij);
                end
                
                % Set vertical positions
                tmp = (0:M-1) * lambda*spacing;
                posv = tmp - mean(tmp);
                tmp = reshape( posv(ones(1,N),:).' , 1 , [] );
                h_array.element_position(3,:) =...
                    reshape( tmp(ones(T,1),:) ,1,[] );
                
                % Set horizontal positions
                tmp = (0:N-1) * lambda*spacing;
                posh = tmp - mean(tmp);
                tmp = reshape( posh(ones(1,M),:) , 1 , [] );
                h_array.element_position(2,:) =...
                    reshape( tmp(ones(T,1),:) ,1,[] );
            end
            
            
        case 'parametric'
            
            % Set inputs if not given
            if ~exist('Ain','var')
                Ain = 1.9;
            end
            if ~exist('Bin','var')
                Bin = 0.1;
            end
            if ~exist('Cin','var')
                Cin = 1;
            end
            if ~exist('Din','var')
                Din = 1.3;
            end
            
            
            if element == 0 || isempty(h_array.no_elements)
                h_array.generate('omni');
                h_array.interpolation_method        = 'linear';
                e = 1;
            else
                e = element;
            end
            
            phi = h_array.azimuth_grid;
            theta = h_array.elevation_grid;
            
            C = cos(theta).^Cin;
            D = exp(-Din * phi.^2);
            
            P = zeros(numel( theta ),numel(phi));
            for a = 1:numel( theta )
                for b = 1:numel(phi)
                    P(a, b) = C(a) * D(b);
                end
            end
            P = Bin + (1-Bin)*P;
            
            h_array.Fa(:,:,e) = Ain * sqrt(P);
            h_array.Fb(:,:,e) = zeros(h_array.no_el, h_array.no_az);
            
        case 'xpol'
            if element ~= 0
                error('Choosing one element is not allowed when generating antennas with more than one element.')
            end
            h_array.generate('omni');
            h_array.copy_element(1,2);
            h_array.rotate_pattern(90,'x',2);
            
        case 'rhcp-dipole'
            if element ~= 0
                error('Choosing one element is not allowed when generating antennas with more than one element.')
            end
            h_array.generate('dipole');
            h_array.generate('dipole',2);
            h_array.rotate_pattern(90,'x',2);
            h_array.coupling = 1/sqrt(2) * [1;-1j];
            
        case 'lhcp-dipole'
            h_array.generate('rhcp-dipole',element);
            h_array.coupling = 1/sqrt(2) * [1;1j];
            
        case 'lhcp-rhcp-dipole'
            h_array.generate('rhcp-dipole',element);
            h_array.coupling = 1/sqrt(2) * [1 1;1j -1j];
            
        case 'ula2'
            if element ~= 0
                error('Choosing one element is not allowed when generating arrays.')
            end
            h_array.generate('omni');
            h_array.no_elements                 = 2;
            h_array.element_position(2,:)       = [-0.05 0.05];
            
        case 'ula4'
            if element ~= 0
                error('Choosing one element is not allowed when generating arrays.')
            end
            h_array.generate('omni');
            h_array.no_elements                 = 4;
            h_array.element_position(2,:)       = -0.15 :0.1: 0.15;
            
        case 'ula8'
            if element ~= 0
                error('Choosing one element is not allowed when generating arrays.')
            end
            h_array.generate('omni');
            h_array.no_elements                 = 8;
            h_array.element_position(2,:)       = -0.35 :0.1: 0.35;
            
    end
    
    h_array.name = array_type;
end
end

