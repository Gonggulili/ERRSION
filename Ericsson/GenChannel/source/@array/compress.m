function ratio = compress( h_array )
%COMPRESS Stores the array in compressed form
%
%   If there are many similar elements in an array, the memory and storage
%   requirements might be high. Therefore, it is possible to compress
%   the array to save storage space. This is done as follows:
%
%   - Patterns are stored in polar-spheric polarization basis
%   - If multiple elements have the same patterns, the pattern is stored
%     only once.
%   - Patterns are stored in single precision.
%   - If there are complex valued patterns with no imaginary part, they are
%     converted to real values.
%
%   It is recommended to call "compress" before saving an array to disk.
%   Decompressing is done automatically when needed.
%
%   The parameter "ratio" is the compression factor in percent.
%
% QuaDRiGa Copyright (C) 2011-2015 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
%
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% Strictly use polar-spheric base
if ~strcmp( h_array.polarization_basis, 'polar-spheric' )
    h_array.change_pol_basis('polar-spheric');
end

if ~h_array.iscompressed
    
    n_el = h_array.no_elements;
    
    % Determine size before compression
    n_val = h_array.no_az * h_array.no_el * n_el;
    orig_size = 2 * n_val;  % Fa and Fb
    if ~isreal( h_array.Fa )
        orig_size = orig_size + n_val;
    end
    if ~isreal( h_array.Fb )
        orig_size = orig_size + n_val;
    end
    
    % Calculate the compressed patterns
    Pind = zeros(2,n_el);
    [ FCa, Pind(1,:) ] = find_equal_patterns( h_array.Fa );
    [ FCb, Pind(2,:) ] = find_equal_patterns( h_array.Fb );
    
    % Determine size after compression
    n_val = prod( prod( size( FCa ) ) );
    new_size = n_val;
    if ~isreal( FCa )
        new_size = 2*n_val;
    end
    n_val = prod( prod( size( FCb ) ) );
    if ~isreal( FCb )
        new_size = new_size + 2*n_val;
    else
        new_size = new_size + n_val;
    end
    
    % Calculate compression factor and set array to single precision.
    switch h_array.precision
        case 'double'
            ratio = 0.5 * new_size / orig_size * 100;
            h_array.precision = 'single';
            
        case 'single'
            ratio = new_size / orig_size * 100;
    end
    
    % Update array object
    h_array.Pind = Pind;
    h_array.PFa = FCa;
    h_array.PFb = FCb;
else
    ratio = 0;
    warning('QuaDRiGa:Array:compress','Array is already compressed.');
end

end



function [ FC , u ] = find_equal_patterns( F )

FC = F(:,:,1);
u = zeros( 1,size(F,3) );
min_min_val = Inf;

unique_el = 0;
while any( u==0 )
    unique_el = unique_el + 1;
    
    % Store the next unique element
    ih = find( u==0,1 );
    FC(:,:,unique_el) = F(:,:,ih);
    u( ih ) = unique_el;
    
    max_val = max( reshape( abs( FC(:,:,unique_el)) , [],1 ) );
    min_val = sqrt( max_val.^2/10.^7 ); % 70 dB dynamic range
    
    if min_min_val > min_val
        min_min_val = min_val;
    end
    
    ii = u == 0;
    ij = ones( 1,sum(ii) ) * unique_el;
    
    Fd = F(:,:,ii) - FC(:,:,ij);
    Fd = reshape( Fd , [] , sum(ii) );
    
    ik = all( abs( Fd ) < min_val , 1 );
    ii = find( ii );
    
    u(1,ii( ik ) ) = unique_el;
end

% Check if array has only real components
if ~isreal( FC )
    if all( imag( FC(:) ) < min_min_val )
        FC = real( FC );
    end
end

% Use single precision
FC = single( FC );

end
