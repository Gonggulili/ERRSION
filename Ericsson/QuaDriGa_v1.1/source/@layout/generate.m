function h_layout = generate( varargin )
%GENERATE Generates predefined layouts
%
%   GENERATE( 'regular' , no_sites , isd , h_array )
%   generates a new multicell layout using a regular grid of BS positions.
%   Each BS has three sectors. The number of sites can be 1, 7, 19 or 37 -
%   resulting in 3, 21, 57 or 111 sectors, respectively. The 'isd'
%   determines the distance between the BS in m. The antenna array
%   'h_array' is for one sector only. It will be rotated to match the
%   sector orientations. The broadside-direction of the provided antenna
%   must be 0 (facing east).
%
%   GENERATE( layout_type , layout_size , no_tx  )
%   generates a new layout with predefined settings.
%
% QuaDRiGa Copyright (C) 2011-2016 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
%
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

supported_types = {'random','regular','regular_base','regular6'};

if nargin > 0
    layout_type = varargin{1};
end

% Check, if the layout type is supported
if nargin == 0 || ~( ischar(layout_type) && any( strcmpi(layout_type,supported_types)) )
    str = 'Layout type not found. Supported types are: ';
    no = numel(supported_types);
    for n = 1:no
        str = [str,supported_types{n}];
        if n<no
            str = [str,', '];
        end
    end
    error(str);
end


switch layout_type
    case 'regular_base'
        % Read the number of sites
        if nargin >= 2
            sites = varargin{2};
            if ~any(sites == [1,7,19,37])
                error('The number of sites must be 1,7, 19 or 37');
            end
        else
            sites = 7;
        end
        
        % Read the ISD
        if nargin >= 3
            isd = varargin{3};
        else
            isd = 500;
        end
        
        h_layout = layout;
        h_layout.no_tx = sites;
        
        % Set the positions of the BS sites
        if sites > 1
            tmp1 = 30+(0:5)*60;
            tmp1 = exp(1j*tmp1*pi/180);
            h_layout.tx_position(1,2:7) = real( tmp1 );
            h_layout.tx_position(2,2:7) = imag( tmp1 );
        end
        if sites > 7
            tmp21 = 2*tmp1;
            h_layout.tx_position(1,9:2:19) = real( tmp21 );
            h_layout.tx_position(2,9:2:19) = imag( tmp21 );
            dist = real(tmp21(1));
            tmp22 = (0:5)*60;
            tmp22 = exp(1j*tmp22*pi/180)*dist;
            h_layout.tx_position(1,8:2:18) = real( tmp22 );
            h_layout.tx_position(2,8:2:18) = imag( tmp22 );
        end
        if sites > 19
            tmp31 = 3*tmp1;
            h_layout.tx_position(1,22:3:37) = real( tmp31 );
            h_layout.tx_position(2,22:3:37) = imag( tmp31 );
            
            tmp32 = real(tmp31(1)) - 1j*0.5;
            tmp32 = tmp32.*exp( 1j*(0:5)*60*pi/180 );
            h_layout.tx_position(1,20:3:35) = real( tmp32 );
            h_layout.tx_position(2,20:3:35) = imag( tmp32 );
            
            tmp33 = real(tmp31(1)) + 1j*0.5;
            tmp33 = tmp33.*exp( 1j*(0:5)*60*pi/180 );
            h_layout.tx_position(1,21:3:36) = real( tmp33 );
            h_layout.tx_position(2,21:3:36) = imag( tmp33 );
        end
        h_layout.tx_position = h_layout.tx_position * isd;
        h_layout.tx_position(3,:) = 25;
        
    case 'regular'
        % Orientations [deg]: 30, 150, -90
        
        h_layout = layout.generate('regular_base',varargin{2:end});
        
        % Read the antenna
        if nargin >= 4
            sector_ant = varargin{4};
            if ~isa(sector_ant,'array')
                error('The tx-antenna must be an array-class object.');
            end
        else
            sector_ant = array;
        end
        
        ant = sector_ant.copy;
        
        % Rotate antenna to match the sector orientation
        no_el = ant.no_elements;
        for n = 1 : no_el
            ant.copy_element( n , n + [1,2]*no_el );
        end
        ant.rotate_pattern( 30 , 'z' , 1:no_el );
        ant.rotate_pattern( 150 , 'z' , no_el+1 : 2*no_el );
        ant.rotate_pattern( -90 , 'z' , 2*no_el+1 : 3*no_el );
        
        h_layout.tx_array = ant;
        
    case 'regular6'
        % Same as regular, but with 6 sectors per BS
        % Orientations [deg]: 0, 60, 120, 180, -120, -60
        
        h_layout = layout.generate('regular_base',varargin{2:end});
        
        % Read the antenna
        if nargin >= 4
            sector_ant = varargin{4};
            if ~isa(sector_ant,'array')
                error('The tx-antenna must be an array-class object.');
            end
        else
            sector_ant = array;
        end
        
        ant = sector_ant.copy;
        
        % Rotate antenna to match the sector orientation
        no_el = ant.no_elements;
        for n = 1 : no_el
            ant.copy_element( n , n + (1:5)*no_el );
        end
        ant.rotate_pattern( 60   , 'z' , no_el+1 : 2*no_el );
        ant.rotate_pattern( 120  , 'z' , 2*no_el+1 : 3*no_el );
        ant.rotate_pattern( 180  , 'z' , 3*no_el+1 : 4*no_el );
        ant.rotate_pattern( -120 , 'z' , 4*no_el+1 : 5*no_el );
        ant.rotate_pattern( -60  , 'z' , 5*no_el+1 : 6*no_el );
        
        h_layout.tx_array = ant;
        
    case 'random'
        
        if nargin >= 2
            layout_size = varargin{2};
        else
            layout_size = 500;
        end
        
        if nargin >= 3
            no_tx = varargin{3};
        else
            no_tx = 1;
        end
        
        
        if ~( isnumeric(layout_size) &&...
                isreal(layout_size) &&...
                all( size(layout_size) == 1 ) &&...
                layout_size > 0 )
            error('??? "layout_size" must be a real scalar > 0');
        end
        
        h_layout = layout;
        h_layout.no_tx = no_tx;
        
        tx_position_new = h_layout.tx_position;
        
        for n = 1:h_layout.no_tx
            a = rand*layout_size * exp( 2*pi*1j * rand );
            b = rand * (40 - 20) + 20;
            
            tx_position_new(1,n) = real(a);
            tx_position_new(2,n) = imag(a);
            tx_position_new(3,n) = b;
            
            h_layout.tx_position = tx_position_new;
        end
        
        h_layout.tx_array(1) = array('dipole');
        for n = 1:h_layout.no_tx
            h_layout.tx_array(n) = h_layout.tx_array(1);
        end
        
        h_layout.rx_array(1) = array('dipole');
        for n = 1:h_layout.no_rx
            h_layout.rx_array(n) = h_layout.rx_array(1);
        end
end

h_layout.name = layout_type;
end
