function change_pol_basis( h_array, new_basis )
%CHANGE_POL_BASE Changes the polarization base of a pattern
%
%   This method can be used to change the polarization basis of an antenna
%   pattern. By default, QuaDRiGa uses the polar-spheric basis. However,
%   the antenna patterns can be given in other bases as well.
%
% QuaDRiGa Copyright (C) 2011-2014 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
%
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% Check if the new basis is valid
if ~array.supported_pol_basis( new_basis )
    str = 'Polarization basis not supported.';
    error('QuaDRiGa:Array:wrongPolarizationBasis',str);
end

% Get the old basis
old_basis   = h_array.polarization_basis;
phi         = h_array.azimuth_grid;
theta       = h_array.elevation_grid';
no_az       = h_array.no_az;
no_el       = h_array.no_el;
no_elements = h_array.no_elements;

% Step 1 - Transform the old basis to 'polar-spheric'.
switch old_basis
    case 'cartesian'
        
        % Transformation matrix
        T = zeros( no_el,no_az,2,3 );
        T(:,:,1,1) =  sin( theta ) * cos( phi );
        T(:,:,1,2) =  sin( theta ) * sin( phi );
        T(:,:,1,3) = -cos( theta ) * ones(1,no_az);
        T(:,:,2,1) =  ones(no_el,1) * -sin( phi );
        T(:,:,2,2) =  ones(no_el,1) * cos( phi );
    
        Fv = zeros( no_el , no_az , no_elements );
        Fh = zeros( no_el , no_az , no_elements );
        
        % Transform Basis
        for i_el = 1 : no_elements
            Fv( :,:,i_el ) = T(:,:,1,1 ) .* h_array.Fa( :,:,i_el ) + ...
                T(:,:,1,2 ) .* h_array.Fb( :,:,i_el ) + ...
                T(:,:,1,3 ) .* h_array.Fc( :,:,i_el );
            
            Fh( :,:,i_el ) = T(:,:,2,1 ) .* h_array.Fa( :,:,i_el ) +...
                T(:,:,2,2 ) .* h_array.Fb( :,:,i_el );
        end
        
        h_array.Fa = Fv;
        h_array.Fb = Fh;
        h_array.Fc = 0;
        h_array.polarization_basis = 'polar-spheric';
        
    case 'az-el'

        % Remove singularity
        ind         = abs( theta ) < 0.001*pi;
        theta(ind)  = 0.001*pi;
        
        sc = sqrt( 1 - cos(theta).^2 * sin(phi).^2 );
        sc = sc ./ ( ones(no_el,1) * cos(phi).^2 + sin(theta).^2 * sin(phi).^2 );

        % Transformation matrix
        T = zeros( no_el,no_az,2,2 );
        T(:,:,1,1) =  ones(no_el,1) * cos( phi );
        T(:,:,2,1) =  -sin( theta ) * sin( phi );
        T(:,:,1,2) =   sin( theta ) * sin( phi );
        T(:,:,2,2) =  ones(no_el,1) * cos( phi );
        T = T .* sc(:,:,[1 1],[1 1]);
              
        Fv = zeros( no_el , no_az , no_elements );
        Fh = zeros( no_el , no_az , no_elements );
        
        % Transform Basis
        for i_el = 1 : no_elements
            Fv( :,:,i_el ) = T(:,:,1,1 ) .* h_array.Fa( :,:,i_el ) +...
                T(:,:,1,2 ) .* h_array.Fb( :,:,i_el );
            
            Fh( :,:,i_el ) = T(:,:,2,1 ) .* h_array.Fa( :,:,i_el ) +...
                T(:,:,2,2 ) .* h_array.Fb( :,:,i_el );
        end
        
        h_array.Fa = Fv;
        h_array.Fb = Fh;
        h_array.Fc = 0;
        h_array.polarization_basis = 'polar-spheric';
        
    case 'el-az'

        % Remove singularity
        ind         = abs( theta ) < 0.001*pi;
        theta(ind)  = 0.001*pi;
        
        sc = sqrt( 1 - cos(theta).^2 * cos(phi).^2 );
        sc = sc ./ ( ones(no_el,1) * sin(phi).^2 + sin(theta).^2 * cos(phi).^2 );

        % Transformation matrix
        T = zeros( no_el,no_az,2,2 );
        T(:,:,1,1) =   sin( theta ) * cos( phi );
        T(:,:,2,1) = -ones(no_el,1) * sin( phi );
        T(:,:,1,2) =  ones(no_el,1) * sin( phi );
        T(:,:,2,2) =   sin( theta ) * cos( phi );
        T = T .* sc(:,:,[1 1],[1 1]);
              
        Fv = zeros( no_el , no_az , no_elements );
        Fh = zeros( no_el , no_az , no_elements );
        
        % Transform Basis
        for i_el = 1 : no_elements
            Fv( :,:,i_el ) = T(:,:,1,1 ) .* h_array.Fa( :,:,i_el ) +...
                T(:,:,1,2 ) .* h_array.Fb( :,:,i_el );
            
            Fh( :,:,i_el ) = T(:,:,2,1 ) .* h_array.Fa( :,:,i_el ) +...
                T(:,:,2,2 ) .* h_array.Fb( :,:,i_el );
        end
        
        h_array.Fa = Fv;
        h_array.Fb = Fh;
        h_array.Fc = 0;
        h_array.polarization_basis = 'polar-spheric';
end


% Step 2 - Transform the antenna to the new basis.
switch new_basis
    case 'cartesian'
        
        % Transformation matrix
        T = zeros( no_el,no_az,3,2 );
        T(:,:,1,1) =  sin( theta ) * cos( phi );
        T(:,:,2,1) =  sin( theta ) * sin( phi );
        T(:,:,3,1) = -cos( theta ) * ones(1,no_az);
        T(:,:,1,2) =  ones(no_el,1) * -sin( phi );
        T(:,:,2,2) =  ones(no_el,1) * cos( phi );
        
        Fx = zeros( no_el , no_az , no_elements );
        Fy = zeros( no_el , no_az , no_elements );
        Fz = zeros( no_el , no_az , no_elements );
        
        % Transform Basis
        for i_el = 1 : no_elements
            Fx( :,:,i_el ) = T(:,:,1,1 ) .* h_array.Fa( :,:,i_el ) +...
                T(:,:,1,2 ) .* h_array.Fb( :,:,i_el );
            
            Fy( :,:,i_el ) = T(:,:,2,1 ) .* h_array.Fa( :,:,i_el ) +...
                T(:,:,2,2 ) .* h_array.Fb( :,:,i_el );
            
            Fz( :,:,i_el ) = T(:,:,3,1 ) .* h_array.Fa( :,:,i_el );
        end
        
        h_array.Fa = Fx;
        h_array.Fb = Fy;
        h_array.Fc = Fz;
        h_array.polarization_basis = 'cartesian';
        
    case 'az-el'
        
        % Remove singularity
        ind         = abs( theta ) < 0.001*pi;
        theta(ind)  = 0.001*pi;
        sc          = sqrt( 1 - cos(theta).^2 * sin(phi).^2 );

        % Transformation matrix
        T = zeros( no_el,no_az,2,2 );
        T(:,:,1,1) =  ones(no_el,1) * cos( phi );
        T(:,:,2,1) =   sin( theta ) * sin( phi );
        T(:,:,1,2) =  -sin( theta ) * sin( phi );
        T(:,:,2,2) =  ones(no_el,1) * cos( phi );
        T = T./sc(:,:,[1 1],[1 1]);
              
        Faz = zeros( no_el , no_az , no_elements );
        Fel = zeros( no_el , no_az , no_elements );
        
        % Transform Basis
        for i_el = 1 : no_elements
            Faz( :,:,i_el ) = T(:,:,1,1 ) .* h_array.Fa( :,:,i_el ) +...
                T(:,:,1,2 ) .* h_array.Fb( :,:,i_el );
            
            Fel( :,:,i_el ) = T(:,:,2,1 ) .* h_array.Fa( :,:,i_el ) +...
                T(:,:,2,2 ) .* h_array.Fb( :,:,i_el );
        end
        
        h_array.Fa = Faz;
        h_array.Fb = Fel;
        h_array.Fc = 0;
        h_array.polarization_basis = 'az-el';
        
    case 'el-az'
        
        % Remove singularity
        ind         = abs( theta ) < 0.001*pi;
        theta(ind)  = 0.001*pi;
        
        sc          = sqrt( 1 - cos(theta).^2 * cos(phi).^2 );

        % Transformation matrix
        T = zeros( no_el,no_az,2,2 );
        T(:,:,1,1) =   sin( theta ) * cos( phi );
        T(:,:,2,1) =   ones(no_el,1) * sin( phi );
        T(:,:,1,2) =  -ones(no_el,1) * sin( phi );
        T(:,:,2,2) =   sin( theta ) * cos( phi );
        T = T./sc(:,:,[1 1],[1 1]);
              
        Faz = zeros( no_el , no_az , no_elements );
        Fel = zeros( no_el , no_az , no_elements );
        
        % Transform Basis
        for i_el = 1 : no_elements
            Faz( :,:,i_el ) = T(:,:,1,1 ) .* h_array.Fa( :,:,i_el ) +...
                T(:,:,1,2 ) .* h_array.Fb( :,:,i_el );
            
            Fel( :,:,i_el ) = T(:,:,2,1 ) .* h_array.Fa( :,:,i_el ) +...
                T(:,:,2,2 ) .* h_array.Fb( :,:,i_el );
        end
        
        h_array.Fa = Faz;
        h_array.Fb = Fel;
        h_array.Fc = 0;
        h_array.polarization_basis = 'el-az';
end

end
