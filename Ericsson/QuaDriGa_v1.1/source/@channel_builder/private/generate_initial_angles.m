function generate_initial_angles( h_cb )
%GENERATE INITIAL_ANGLES Generate angular parameters (private)
%
% In three dimensional space, the signals originating from a cluster would
% arrive at the receiver under a certain angle in azimuth and elevation
% direction. The same hold for the angles of departure.
%
% Note: All generated angles are in radians!
%
% QuaDRiGa Copyright (C) 2011-2014 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
%
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

use_ground_reflection = h_cb.par.simpar.use_ground_reflection;

if use_ground_reflection
    L = h_cb.NumClusters - 1;  % no. taps
else
    L = h_cb.NumClusters;
end
N = h_cb.par.no_positions;                  % no. positions

oL = ones(1,L);
oN = ones(1,N);

angles = h_cb.par.get_angles.'*pi/180;      % Angles between Tx and Rx

if use_ground_reflection
    pow = h_cb.pow(:,[1,3:end]);
    pow(:,1) = pow(:,1) + h_cb.pow(:,2);
else
    pow = h_cb.pow;                         % Path powers
end

spreads = [ h_cb.par.asD ; h_cb.par.asA ; ...
    h_cb.par.esD ; h_cb.par.esA ].' * pi/180;

drifting_precision = h_cb.par.simpar.drifting_precision;

% Azimuth spreads
if L > 1
    for i_spread = 1:4
        
        as    = spreads(:,i_spread);
        sigma = zeros( N,L );
        aso   = zeros( N,1 );
        idx   = true(N,1);
        
        cnt = 1;
        while any( idx ) && cnt < 20
            
            sigma( idx,: ) = as(idx,oL);
            sigma( idx,: ) = sigma( idx,: ) .* randn(sum( idx ),L);
            sigma( idx,1 ) = 0;
            
            if i_spread > 2 % Elevation angles
                % For the elevation angles, we have to include the LOS
                % elevation.
                sigma = sigma + angles(:,i_spread*oL);
            end
            
            tmp = sigma( idx,: );
            for n = 1:5
                
                % Map angels to the allowed range if they are outside
                tmp = mod( real(tmp) + pi, 2*pi) - pi;
                
                if i_spread > 2 % Elevation angles
                    tmp( tmp >  pi/2 ) =  pi - tmp( tmp >  pi/2 );
                    tmp( tmp < -pi/2 ) = -pi - tmp( tmp < -pi/2 );
                end
                
                aso( idx,: ) = qf.calc_angular_spreads( tmp, pow( idx,: ) );
                scale        = as( idx,: )./aso( idx,: );
                
                % Limit the scaling factor to 100 in order to prevent
                % numeric errors
                scale( scale > 100 ) = 100;
                scale( isnan(scale) ) = 1; % Spread = 0
                
                % Determine the angles that are outside the possible range
                if i_spread > 2 % Elevation angles
                    
                    tmp      = scale(:,oL) .* (tmp - angles(idx,i_spread*oL)) +...
                        angles(idx,i_spread*oL);
                    
                    out      = tmp > pi/2;
                    tmp(out) = randn( sum(out(:)) , 1 ) * pi/4 + pi/2;
                    
                    out      = tmp < -pi/2;
                    tmp(out) = randn( sum(out(:)) , 1 ) * pi/4 - pi/2;
                    
                else % Azimuth angles
                    
                    tmp      = scale(:,oL) .* tmp;
                    out      = tmp > pi | tmp < -pi;
                    tmp(out) = randn( sum(out(:)) , 1 ) * pi/2 + pi;
                    
                end
            end
            
            % Map angels to the allowed range if they are outside
            tmp = mod( real(tmp) + pi, 2*pi) - pi;
            if i_spread > 2 % Elevation angles
                tmp( tmp >  pi/2 ) =  pi - tmp( tmp >  pi/2 );
                tmp( tmp < -pi/2 ) = -pi - tmp( tmp < -pi/2 );
            end
            sigma( idx,: ) = tmp;
            
            aso( idx,: )   = qf.calc_angular_spreads( sigma( idx,: ), pow( idx,: ) );
            
            if cnt == 1
                sigma_old = sigma;
                aso_old   = aso;
            else
                id_better = aso > aso_old;
                id_worse  = aso < aso_old;
                
                sigma_old( id_better,: ) = sigma( id_better,: );
                sigma( id_worse,: )      = sigma_old( id_worse,: );
                
                aso_old( id_better ) = aso( id_better );
                aso( id_worse )      = aso_old( id_worse );
            end
            
            idx = aso < 0.95*as | aso > 1.05*as;
            cnt = cnt + 1;
        end
        
        if i_spread <= 2 % Azimuth angles
            sigma = sigma + angles(:,i_spread*oL);
            sigma = mod( sigma + pi, 2*pi) - pi;
        end
        
        switch i_spread
            case 1
                h_cb.AoD = sigma;
            case 2
                h_cb.AoA = sigma;
            case 3
                h_cb.EoD = sigma;
            case 4
                h_cb.EoA = sigma;
        end
    end
else
    h_cb.AoD = angles(:,1);
    h_cb.AoA = angles(:,2);
    h_cb.EoD = angles(:,3);
    h_cb.EoA = angles(:,4);
end

if use_ground_reflection
    
    % Calculate the angles for the ground reflection
    tmp = h_cb.par.positions(3,:);
    h_cb.par.positions(3,:) = -tmp;
    angles_gr = h_cb.par.get_angles.'*pi/180;      % Angles for the ground reflection
    
    h_cb.par.positions(3,:) = tmp;
    
    % Add the ground reflection angles to the list of angles
    AoD = h_cb.AoD( :,[1,1,2:end] );
    AoA = h_cb.AoA( :,[1,1,2:end] );
    EoD = [ h_cb.EoD(:,1), angles_gr(:,3), h_cb.EoD(:,2:end) ];
    EoA = [ h_cb.EoA(:,1), -angles_gr(:,4), h_cb.EoA(:,2:end) ];
    pow = h_cb.pow;
    fixed = [ true(N,2), false(N,L-1) ];
    
    % Correct the NLOS angles to optain the given angular spread as good as possible
    h_cb.AoD = update_angles( spreads(:,1), AoD, pow, fixed, [-pi,pi] );
    h_cb.AoA = update_angles( spreads(:,2), AoA, pow, fixed, [-pi,pi] );
    h_cb.EoD = update_angles( spreads(:,3), EoD, pow, fixed, [-pi,pi]/2 );
    h_cb.EoA = update_angles( spreads(:,4), EoA, pow, fixed, [-pi,pi]/2 );
    
    L = L+1;
    oL = ones(1,L);
end

if L > 1
    if drifting_precision == 4
        % In this mode, the departure angles depend on the arrival angles.
        % Hence, we calculate the corresponding departure angles.
        
        % Calculate the unit length vector pointing from the Rx to the LBS.
        [ ahat_x, ahat_y, ahat_z ] = sph2cart( h_cb.AoA, h_cb.EoA, 1);
        ahat = permute( cat( 3, ahat_x, ahat_y, ahat_z ) , [3,1,2] );
        
        r = cat( 2,h_cb.par.rx_track.initial_position ) - h_cb.par.tx_position(:,oN);
        norm_r = sqrt(sum(r.^2)).';
        
        dist = h_cb.taus.*simulation_parameters.speed_of_light + norm_r(:,oL);
        
        [ b, norm_b ]  = solve_cos_theorem( ahat , r , dist  );
        
        h_cb.AoD = permute( atan2(b(2,:,:),b(1,:,:) ) , [2,3,1]);
        h_cb.EoD = permute( asin( b(3,:,:)./norm_b  ) , [2,3,1]);
        
        
    elseif drifting_precision == 5
        % This implements the multi-bounce model.
        % Those values are constant and need to be calculated only once
        
        [ ahat_x, ahat_y, ahat_z ] = sph2cart(h_cb.AoA, h_cb.EoA, 1);
        ahat = permute( cat( 3, ahat_x, ahat_y, ahat_z ) , [3,1,2] );
        ahat = reshape( ahat , 3 , N*L );
        
        [ bhat_x, bhat_y, bhat_z ] = sph2cart(h_cb.AoD, h_cb.EoD, 1);
        bhat = permute( cat( 3, bhat_x, bhat_y, bhat_z ) , [3,1,2] );
        bhat = reshape( bhat , 3 , N*L );
        
        r = cat( 2,h_cb.par.rx_track.initial_position ) - h_cb.par.tx_position(:,oN);
        norm_r = sqrt(sum(r.^2)).';
        r = reshape( r(:,:,oL) , 3 , N*L );
        
        dist = h_cb.taus.*simulation_parameters.speed_of_light + norm_r(:,oL);
        dist = permute( dist , [3,1,2] );
        dist = reshape( dist , 1 , N*L );
        
        % Test if the optimization problem can be solved and update the angles
        % (single-bounce) that violate the assumptions.
        [ ~, ~, ~, valid ] = solve_multi_bounce_opti( ahat, bhat, r, dist, 2 );
        invalid = reshape( ~valid, N,L ) ;
        invalid( :,1 ) = false;  % LOS is always valid
        
        if use_ground_reflection
            invalid( :,2 ) = false;  % Ground Reflection is always valid
        end
        
        % Use the single-bounce model for points where the problem has no
        % solution.
        [ b, norm_b ] = solve_cos_theorem( ahat(:,invalid(:)),...
            r(:,invalid(:)), dist(invalid(:)));
        h_cb.AoD(invalid) = atan2(b(2,:,:),b(1,:,:) );
        h_cb.EoD(invalid) = asin( b(3,:,:)./norm_b  );
        
        % Here we rescale the angles in order to maintain the angular spreads.
        rescale = find( sum( invalid ,2 ) > 0 );
        if any( rescale )
            fixed = invalid( rescale,: );
            if use_ground_reflection
                fixed(:,[1,2]) = true;
            else
                fixed(:,1) = true;
            end
            
            % AoD rescaling
            AoD = update_angles( spreads(rescale,1), h_cb.AoD(rescale,:), pow(rescale,:), fixed, [-pi,pi] );
            h_cb.AoD( rescale,: ) = AoD;
            
            % EoD rescaling
            EoD = update_angles( spreads(rescale,3), h_cb.EoD(rescale,:), pow(rescale,:), fixed, [-pi/2,pi/2] );
            h_cb.EoD( rescale,: ) = EoD;
        end
    end
end

end
