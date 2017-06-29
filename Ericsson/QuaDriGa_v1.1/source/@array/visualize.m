function h_figures = visualize( h_array, element_indices )
%VISUALIZE Create a plot showing the element configurations
%
%   VISUALIZE plots the 3D beam patterns of all antenna elements within the
%             antenna array.
%   VISUALIZE(element) only plots the beam patterns of the antenna elements
%             specified by the given element indices contained in 'element'.
%
% QuaDRiGa Copyright (C) 2011-2014 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
%
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

if nargin == 1
    element_indices = 1:h_array.no_elements;
end

if ~(any(size(element_indices) == 1) && isnumeric(element_indices) ...
        && isreal(element_indices) && all(mod(element_indices, 1) == 0) && all(element_indices > 0))
    error('??? "element" must be integer and > 0')
elseif any(element_indices > h_array.no_elements)
    error('??? "element" exceed "no_elements"')
end

h_array.uncompress;

az = double(h_array.azimuth_grid);
if numel(az) == 361
    indt = 1:5:361;
else
    indt = 1:numel(az);
end
az = az(indt);

elev = double(h_array.elevation_grid');
if numel(elev) == 181
    indp = 1:5:181;
else
    indp = 1:numel(elev);
end
elev = elev(indp);

[az_grid, elev_grid] = meshgrid(az, elev);
[Xi, Yi, Zi] = sph2cart(az_grid, elev_grid, 1);

min_value = -20;
scaling = 1;

h_figures = zeros(1, numel(element_indices));

for i_element_indices = 1:numel(element_indices)
    
    % Read the array elements
    Fa = h_array.Fa(indp, indt, element_indices(i_element_indices));
    Fb = h_array.Fb(indp, indt, element_indices(i_element_indices));
    if numel( h_array.Fc ) == 1
        Fc = zeros( size(Fa) );
    else
        Fc = h_array.Fc(indp, indt, element_indices(i_element_indices));
    end
    
    % calculate radiation power pattern
    P = abs(Fa).^2 + abs(Fb).^2 + abs(Fc).^2;
    
    % Normalize by max value
    P_max = max( P(:) );
    P = P ./ P_max;
    
    % Calculate directivity
    tmp         = cos(elev_grid(:));
    directivity_lin    = sum(tmp) ./ sum(P(:).*tmp);
    directivity_dbi    = 10*log10(directivity_lin);
    
    % Normalize patterns by directivity
    tmp = sqrt(directivity_lin./P_max);
    Fa = Fa .* tmp;
    Fb = Fb .* tmp;
    Fc = Fc .* tmp;
    
    h_figures(i_element_indices) = figure('Position', [50, 400, 1200, 500],...
        'Name', [h_array.name ' element ', num2str(element_indices(i_element_indices))]);
    
    axes('position',[0 0 0.92 0.9]);%,'Visible','Off');
    axis off
    
    %title(['Array Antenna Element ', num2str(element_indices(i_element_indices))] );
    
    switch h_array.polarization_basis
        case 'cartesian' 
            text(0.195,0.96,'D^{[x]}(\theta, \phi)','HorizontalAlignment','center');
            text(0.515,0.96,'D^{[y]}(\theta, \phi)','HorizontalAlignment','center');
            text(0.84 ,0.96,'D^{[z]}(\theta, \phi)','HorizontalAlignment','center');
        case 'polar-spheric' 
            text(0.35,1,'D^{[\theta]}(\theta, \phi)');
            %text(0.9 ,1,'D^{[\phi]}(\theta, \phi)');
        case 'az-el' 
            text(0.02,1,'D^{[Az]}(\theta, \phi)');
            text(0.9 ,1,'D^{[El]}(\theta, \phi)');
        case 'el-az' 
            text(0.02,1,'D^{[\alpha]}(\theta, \phi)');
            text(0.9 ,1,'D^{[\epsilon]}(\theta, \phi)');
    end
    
    if strcmp( h_array.polarization_basis ,'cartesian' )
        no_plots = 3;
    else
        no_plots = 1;
    end
    
    for m = 1:no_plots
        if no_plots == 3
            axes('position',[-0.25+0.3*m 0.12 0.25 0.7]);
            Po = 10*log10(abs(Fa).^2);
        else
            axes('position',[-0.1+0.45*m 0.12 0.38 0.82]);
            Po = 10*log10(abs(Fa).^2);
        end
        
        switch m
            case 1
                Po = 10*log10( abs(Fa).^2 );
            case 2
                Po = 10*log10( abs(Fb).^2 );
            case 3
                Po = 10*log10( abs(Fc).^2 );
        end
        Po = double( Po );
        
        P = Po;
        P(P < min_value) = min_value;
        P = (P - min_value) ./ (directivity_dbi - min_value) .* scaling;
        
        X = P .* Xi;
        Y = P .* Yi;
        Z = P .* Zi;
        
        surf(X, Y, Z, Po)
        
        axis equal
        axis(scaling.*[-1 1 -1 1 -1 1]);
        caxis([min_value, directivity_dbi]);
        set(gca, 'xtick', (-1:1).*scaling/2);
        set(gca, 'ytick', (-1:1).*scaling/2);
        set(gca, 'ztick', (-1:1).*scaling/2);
        xlabel('x')
        ylabel('y')
        zlabel('z')
       
        view(45, 33)
    end
    
    axes('position', [0.08 0.08 0.08 0.08],'Visible','Off');
    caxis([min_value, directivity_dbi]);
    han = colorbar('EastOutside', 'XTick', min_value:3:floor(directivity_dbi));
    set(han, 'position', [0.72 0.06 0.02 0.90])
    zlab = get(han, 'ylabel');
    set(zlab, 'String', 'Partial Directivity in dBi');
    
end


end
