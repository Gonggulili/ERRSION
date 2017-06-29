function aso = init_parameters( h_cb, force )
%INIT_PARAMETERS Generates the initial parameters
%
%  This function creates the initial parameters for the channel builder. If
%  the parameters are already initialized, no new update is performed. The
%  optional parameter 'force' can be used to enforce an update, even if the
%  parameters are already given.
%
% Output:
%   aso     The obtained angular spread in [deg]
%
% QuaDRiGa Copyright (C) 2011-2016 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
%
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

if ~exist( 'force' , 'var' )
    force = false;
end

% Set the number of clusters depending on the frequency
if h_cb.par.scenpar.NumClusters_gamma == 0
    n_clusters = h_cb.par.scenpar.NumClusters;
else
    n_clusters = ceil( h_cb.par.scenpar.NumClusters + ...
        h_cb.par.scenpar.NumClusters_gamma .* ...
        log10(h_cb.par.simpar.center_frequency/1e9) );
    n_clusters( n_clusters<1 ) = 1;
end

if h_cb.par.simpar.use_ground_reflection
    n_clusters = n_clusters + 1;
end

h_cb.NumClusters = n_clusters;

n_subpaths      = 20;
n_paths         = n_clusters*n_subpaths;
n_mobiles       = h_cb.par.no_positions;

if isempty(h_cb.taus) || force
    generate_initial_paths( h_cb );
end

if isempty(h_cb.AoD) || force
    switch h_cb.par.simpar.use_angular_mapping
        case 1
            generate_initial_angles_winner( h_cb );
        case 2
            generate_initial_angles( h_cb );
    end
end

if nargout > 0
    N = size( h_cb.pow,1 );
    aso = zeros( 4,N );
    aso( 1 , : ) = calc_angular_spreads( h_cb.AoD,  h_cb.pow )*180/pi;
    aso( 2 , : ) = calc_angular_spreads( h_cb.AoA,  h_cb.pow )*180/pi;
    aso( 3 , : ) = calc_angular_spreads( h_cb.EoD,  h_cb.pow )*180/pi;
    aso( 4 , : ) = calc_angular_spreads( h_cb.EoA,  h_cb.pow )*180/pi;
end

if isempty(h_cb.xpr) || force
    generate_initial_xpr( h_cb );
end

% Random initial phases
if isempty(h_cb.pin) || force
    pin = ( rand(n_mobiles,n_paths)-0.5 )*2*pi;
    pin = reshape( pin, n_mobiles, n_clusters, n_subpaths );
    
    % LOS path has zero-phase
    pin(:,1,:) = 0;
    
    if h_cb.par.simpar.use_ground_reflection
        % Ground reflection has phase offset of 180 deg
        pin(:,2,:) = pi;
    end
    
    h_cb.pin = reshape(pin,n_mobiles,[]);
end

if isempty(h_cb.subpath_coupling) || force
    [~,h_cb.subpath_coupling] =...
        sort(rand(n_subpaths,4,n_clusters,n_mobiles));
end

end
