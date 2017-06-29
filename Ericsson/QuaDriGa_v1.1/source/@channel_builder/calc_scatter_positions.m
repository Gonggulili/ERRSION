function [ phi_d_lm, theta_d_lm , lbs_pos , fbs_pos ] =...
    calc_scatter_positions( h_cb, i_mobile )
%CALC_SCATTER_POSITIONS Calculates the positions of the scatterers
%
%   This function calculates the positions of the scatterers and
%   initializes the drifting module. The output variables are the NLOS Tx
%   angles for the pre-computation of the Tx array response. 
%
% Input:
%   i_mobile    The index of the mobile terminal.
%
% Output:
%   phi_d_lm    The departure azimuth angles.
%   theta_d_lm  The departure elevation angles.
%   lbs_pos     The position of the last-bounce scatterer.
%   fbs_pos     The position of the first-bounce scatterer.
%
% QuaDRiGa Copyright (C) 2011-2014 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
%
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% Read some common variables
n_path              = h_cb.NumClusters;
n_subpath           = 20;
o_path              = ones(1,n_path);
o_subpath           = ones(1,n_subpath);
drifting_precision  = h_cb.par.simpar.drifting_precision;

% The distance vector from the Tx to the initial position
r = h_cb.par.rx_track(i_mobile).initial_position - h_cb.par.tx_position;
norm_r = sqrt(sum(r.^2)).';

% Get the total path length of the NLOS component
dist = h_cb.taus(i_mobile, :).*simulation_parameters.speed_of_light + norm_r(1,o_path);
dist = dist(:,:,o_subpath);

% We get the initial angles with random coupling
[ phi_d_lm, theta_d_lm, phi_a_lm, theta_a_lm ] =...
    h_cb.get_subpath_angles(i_mobile);

% Get the direction of the last bounce scatterer (LBS) seen from the
% receivers initial position.
[ ahat_lm_x, ahat_lm_y, ahat_lm_z ] = sph2cart(phi_a_lm, theta_a_lm, 1);
ahat_lm = [ ahat_lm_x; ahat_lm_y; ahat_lm_z ];

if drifting_precision == 5
    % Multi-Bounce model.
    
    % Get the direction of the first bounce scatterer (FBS) seen from the
    % transmitter center position.
    [ bhat_lm_x, bhat_lm_y, bhat_lm_z ] = sph2cart(phi_d_lm, theta_d_lm, 1);
    bhat_lm = [ bhat_lm_x; bhat_lm_y; bhat_lm_z ];
    
    [ norm_a_lm, norm_b_lm,~,valid ] =...
        solve_multi_bounce_opti( ahat_lm, bhat_lm, r, dist, 2 );

    % LOS
    norm_a_lm(1,1,:) = 0.5*norm_r;
    norm_b_lm(1,1,:) = 0.5*norm_r;
    valid(1,1,:) = true;
    
    % For the invalid paths, use single-bounce-model
    iv = ~valid;
    if any( iv(:) )
        ahat_iv = [ ahat_lm_x(iv) , ahat_lm_y(iv) , ahat_lm_z(iv) ].';
        [ b, norm_b_lm(iv), norm_a_lm(iv) ]  =...
            solve_cos_theorem( ahat_iv , r , dist(iv).'  );
        
        tmp = norm_b_lm(iv).';
        bhat_lm_x(iv) = b(1,:)./tmp;
        bhat_lm_y(iv) = b(2,:)./tmp;
        bhat_lm_z(iv) = b(3,:)./tmp;
    end
    
    % Calculate the FBS position (relative to initial Rx-pos)
    fbs_pos = zeros( 3, n_path, n_subpath  );
    fbs_pos(1,:,:) = norm_b_lm .* bhat_lm_x - r(1);
    fbs_pos(2,:,:) = norm_b_lm .* bhat_lm_y - r(2);
    fbs_pos(3,:,:) = norm_b_lm .* bhat_lm_z - r(3);
    
else
    % Single-Bounce model
    [ ~, norm_b_lm, norm_a_lm ] = solve_cos_theorem( ahat_lm , r , dist  );
    
    % LOS
    norm_a_lm(1,1,:) = 0.5*norm_r;
    norm_b_lm(1,1,:) = 0.5*norm_r;
end

% Calculate the LBS position
lbs_pos = zeros( 3, n_path, n_subpath  );
lbs_pos(1,:,:) = norm_a_lm .* ahat_lm_x;
lbs_pos(2,:,:) = norm_a_lm .* ahat_lm_y;
lbs_pos(3,:,:) = norm_a_lm .* ahat_lm_z;

if drifting_precision < 5
   fbs_pos = lbs_pos; 
end

% Here, we calculate the vector pointing from each Tx element to the
% initial Rx position. This vector is stored in "e_tx" and update the
% departure angles.

if drifting_precision >= 4
    e_tx       = h_cb.par.tx_array.element_position;
    n_tx       = h_cb.par.tx_array.no_elements;
    o_tx       = ones(1,n_tx);
    
    r_t        = r(:,o_tx) - e_tx;
    r_t        = permute( r_t,[1,3,4,5,2] );
    
    b_tlm      = r_t(:,o_path,o_subpath,1,:) + fbs_pos(:,:,:,1,o_tx) ;
    norm_bc    = sqrt( sum( b_tlm.^2 , 1 ) );
    
    phi_d_lm   = atan2(b_tlm(2,:,:,:,:), b_tlm(1,:,:,:,:));
    theta_d_lm = asin(b_tlm(3,:,:,:,:)./norm_bc);
    
    if drifting_precision == 5
        norm_c_lm = sqrt(sum((lbs_pos - fbs_pos).^2));
        norm_bc   = norm_bc + norm_c_lm( 1,:,:,1,o_tx );
    end
    
    % Initialize the drifting function
    h_cb.update_drifting( 1, 1, i_mobile, lbs_pos, r_t, norm_bc, r );
else
    h_cb.update_drifting( 1, 1, i_mobile, lbs_pos, r, norm_b_lm, r );
end

% Calculate the scatterer positions in global coordinates.
% The positions are relative to the initial Rx-position. This is corrected
% here.
if nargout > 2
    rx_pos = h_cb.par.rx_track( i_mobile ).initial_position;
    lbs_pos = lbs_pos + rx_pos(:,o_path,o_subpath );
    fbs_pos = fbs_pos + rx_pos(:,o_path,o_subpath );
end

end
