function h_channel = get_channels( h_cb, vb_dots )
%GET_CHANNELS Generates channel coefficients
%
%   h_channel = GET_CHANNELS generates the channel coefficients. This is
%   the main function of the channel builder.
%
% QuaDRiGa Copyright (C) 2011-2016 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
%
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

verbose = h_cb(1).par.simpar.show_progress_bars;
if verbose && nargin == 1
    fprintf('Channels     [');
    vb_dots = 50;
    tStart = clock;
end
m0=0;

if numel(h_cb) > 1
    
    vb_dots = zeros( 1,numel(h_cb) );
    if verbose
        for i_cb = 1:numel(h_cb)
            vb_dots(i_cb) = ...
                h_cb(i_cb).par.no_positions;
        end
        vb_dots = qf.init_progress_dots(vb_dots);
    end
    
    h_channel = channel.empty;
    for i_cb = 1:numel(h_cb)
        tmp = h_cb(i_cb).get_channels(vb_dots(i_cb));
        h_channel = [ h_channel , tmp ];
    end
    
else
    % These variables are often needed. Pre-computing them saves a lot of time
    use_polarization_rotation = h_cb.par.simpar.use_polarization_rotation;
    use_random_initial_phase = h_cb.par.simpar.use_random_initial_phase;
    use_ground_reflection = h_cb.par.simpar.use_ground_reflection;
    map_valid = h_cb.par.map_valid;
    wave_no = 2*pi/h_cb.par.simpar.wavelength;
    initial_rx_position = h_cb.par.positions;
    drifting_precision = h_cb.par.simpar.drifting_precision;
    
    % Access to class-properties is time consuming.
    % The array interpolation is the most time intense operation in the
    % channel builder. We save some computing time by reading the arrays
    % here and passing them as variables to the interpolate function later
    % on.
    tx_elevation_grid   = h_cb.par.tx_array.elevation_grid;
    tx_azimuth_grid     = h_cb.par.tx_array.azimuth_grid;
    tx_patV             = h_cb.par.tx_array.Fa;
    tx_patH             = h_cb.par.tx_array.Fb;
    tx_element_pos      = h_cb.par.tx_array.element_position;
    
    % Generate the initial parameters here
    % If the parameters are already given in the channel builder object,
    % we use the given parameters.
    h_cb.init_parameters;

    % Set initial parameters
    n_clusters      = h_cb.NumClusters;
    n_subpaths      = 20;
    n_paths         = n_clusters*n_subpaths;
    n_mobiles       = h_cb.par.no_positions;
    n_tx            = h_cb.par.tx_array.no_elements;
    
    o_clusters      = ones(1,n_clusters);
    o_subpaths      = ones(1,n_subpaths);
    o_tx            = ones(1,n_tx);
    
    if drifting_precision == 0 || map_valid == 0
        % Get the path loss for each rx position.
        % When drifting is enabled, the path loss is also drifting und will
        % be calculated separately.
        [ path_loss , scale_sf ] = h_cb.par.get_pl;
        rx_power = -path_loss.' + 10*log10( h_cb.par.sf ) .* scale_sf;
        rx_power = sqrt( 10.^( 0.1 * rx_power ) );
    end
    
    % The loop for each user position
    h_channel = channel.empty(n_mobiles,0);
    for i_mobile = 1 : n_mobiles
        if verbose; m1=ceil(i_mobile/n_mobiles*vb_dots); if m1>m0;
                for m2=1:m1-m0; fprintf('o'); end; m0=m1; end;
        end;
        
        % It is possible for all users to share the same track object. In
        % this case, the initial position will be the same for all tracks.
        % However, in the parameter_set object, different initial positions
        % can be specified. The following statement compares both values
        % and corrects the initial position in the tack object, if needed.
        
        if any( h_cb.par.rx_track(i_mobile).initial_position ~= ...
                initial_rx_position(:,i_mobile) )
            h_cb.par.rx_track(i_mobile).initial_position =...
                initial_rx_position(:,i_mobile);
        end
        
        % Read some commonly needed variables in order to save time.
        n_rx        = h_cb.par.rx_array(i_mobile).no_elements;
        n_links     = n_rx*n_tx;
        n_snapshots = h_cb.par.rx_track(i_mobile).no_snapshots;
        initial_pos = h_cb.par.rx_track(i_mobile).segment_index( ...
            min( [h_cb.par.rx_track(i_mobile).no_segments,2] ));
        
        % Access to class-properties is time consuming.
        % The array interpolation is the most time intense operation in the
        % channel builder. We save some computing time by reading the arrays
        % here and passing them as variables to the interpolate function.
        % Since very often, mobiles have identical antennas, we only read
        % the data again, if it changes.
        if i_mobile == 1 || isequal( h_cb.par.rx_array(i_mobile) , h_cb.par.rx_array(i_mobile-1) )
            rx_elevation_grid   = h_cb.par.rx_array(i_mobile).elevation_grid;
            rx_azimuth_grid     = h_cb.par.rx_array(i_mobile).azimuth_grid;
            rx_patV             = h_cb.par.rx_array(i_mobile).Fa;
            rx_patH             = h_cb.par.rx_array(i_mobile).Fb;
            rx_element_pos      = h_cb.par.rx_array(i_mobile).element_position;
        end
        
        % Extract the random initial phases
        pin = h_cb.pin(i_mobile,:);
        
        % We need the directions. If they are not provided, compute them here.
        if isempty( h_cb.par.rx_track(i_mobile).ground_direction )
            h_cb.par.rx_track(i_mobile).compute_directions;
        end
        
        switch drifting_precision
            case 0  % Drifting based on rotating phasors
                
                % If we don't use drifting and have a linear track, then the
                % Doppler component is only dependent on the rotating phases of the
                % taps. So, we don't recalculate the antenna response for each
                % snapshot.
                
                % Get the angles of the 20 subpaths and perform random coupling.
                [ aod,eod,aoa,eoa,delay ] =...
                    h_cb.get_subpath_angles(i_mobile);
                
                % Calculate the distance-dependent phases
                lambda  = h_cb.par.simpar.wavelength;
                r       = h_cb.par.rx_track(i_mobile).initial_position - h_cb.par.tx_position;
                norm_r  = sqrt(sum(r.^2)).';
                d_lms   = norm_r + h_cb.par.simpar.speed_of_light * delay;
                phase   = 2*pi/lambda * mod(d_lms, lambda);
                phase   = phase( :,:,o_subpaths );
                
                % Doppler component
                % Without drifting, the Doppler component is calculated by
                % plane wave approximation using the distance from the initial
                % position.
                tmp = h_cb.par.rx_track(i_mobile).positions;
                dist = sqrt( sum([ tmp(1,:) - tmp(1,1) ; ...
                    tmp(2,:) - tmp(2,1)   ; ...
                    tmp(3,:) - tmp(3,1)   ].^2 ) );
                
                % Generate the tx-array channel coefficients for each user position.
                % Since the receiver is mobile, we have to adjust for the movement direction
                % inside the loop.
                
                [Vt,Ht,Pt] = h_cb.par.tx_array.interpolate( ...
                    aod, eod, 1:n_tx, tx_azimuth_grid, tx_elevation_grid, ...
                    tx_patV, tx_patH, tx_element_pos );
                Pt = reshape( Pt, 1 ,n_paths, n_tx );
                no_snap_process = 1;
                
            case { 1,2,3 }  % Tx uses planar waves
                
                [ aod, eod ] = h_cb.calc_scatter_positions( i_mobile );
                
                % Pre-calculate the Tx-Array-Pattern response for the first snapshot
                [ Vt, Ht, Pt ] = h_cb.par.tx_array.interpolate(...
                    aod, eod, 1:n_tx, tx_azimuth_grid, tx_elevation_grid,...
                    tx_patV, tx_patH, tx_element_pos );
                Pt  = reshape( Pt, 1, n_paths, n_tx );
                
                if drifting_precision == 1
                    delay = zeros( n_snapshots, n_clusters );
                else
                    delay = zeros( n_snapshots, n_clusters, n_rx );
                end
                no_snap_process = n_snapshots;
                
            case {4,5}  % Tx uses spherical waves
                
                [ aod, eod ] = h_cb.calc_scatter_positions( i_mobile );
                
                % Pre-calculate the Tx-Array-Pattern response for each
                % element separately.
                Vt  = zeros(1,n_clusters,n_subpaths,n_tx);
                Ht  = zeros(1,n_clusters,n_subpaths,n_tx);
               
                for i_tx = 1:n_tx
                    [ Vt(1,:,:,i_tx), Ht(1,:,:,i_tx)  ] =...
                        h_cb.par.tx_array.interpolate( ...
                        aod(:,:,:,:,i_tx), eod(:,:,:,:,i_tx), i_tx, ...
                        tx_azimuth_grid, tx_elevation_grid, tx_patV,...
                        tx_patH, tx_element_pos );
                end
                
                delay = zeros( n_snapshots, n_clusters, n_rx, n_tx );
                no_snap_process = n_snapshots;
        end
        
        % Travel directions
        gdir = h_cb.par.rx_track(i_mobile).ground_direction;
        hdir = h_cb.par.rx_track(i_mobile).height_direction;
        
        % XPR-Values and array indices
        if use_polarization_rotation > 0
            xprmat = zeros(4,n_paths);                         % Initialization of the pol. rotation
        end
        if use_polarization_rotation == 0 || use_polarization_rotation == 3
            % This is the polarization coupling from WINNER. The
            % polarization is initialized by random phases which are scaled
            % by the XPR. Polarization drifting is not supported.
            
            % XPR-Values
            xpr = 10.^( h_cb.xpr(i_mobile,:,:)/10 );
            xpr = sqrt( reshape( 1./xpr ,1,n_paths ) );
            
            % Random initial phases
            xprmat_winner = exp( 1j*h_cb.random_pol(:,:,i_mobile));
            
            % Global XPR-Matrix
            xprmat_winner = xprmat_winner .* [ones(1,n_paths);xpr;xpr;ones(1,n_paths)];
            
            % Identity Matrix for the LOS-Path
            xprmat_winner( 1,1:n_clusters:n_paths ) = 1;
            xprmat_winner( 2,1:n_clusters:n_paths ) = 0;
            xprmat_winner( 3,1:n_clusters:n_paths ) = 0;
            xprmat_winner( 4,1:n_clusters:n_paths ) = 1;
            
            % Set rotation angles to 0
            gamma = zeros( 1,n_clusters,n_subpaths );
            
        else % use_polarization_rotation == 1 || use_polarization_rotation == 2
            
            % Conversion of the XPR into one rotation angle
            gamma = acot( sqrt( 10.^(0.1*h_cb.xpr(i_mobile,:,:)) ) );
            gamma = reshape( gamma,1,n_clusters,n_subpaths );             % Rotation angles
            
            if use_polarization_rotation == 2
                kappa = exp(1j*h_cb.kappa(i_mobile,:,:));  % Optional HV phase offset
                kappa = reshape(kappa,1,n_paths);
            end
        end
        
        % Placeholder for the coefficient calculation
        cn    = zeros( n_links , n_clusters , n_snapshots );
        
        % Placeholder for the radiated power
        ppat  = zeros( n_links , n_clusters , n_snapshots );
        
        % Do for each snapshot
        for i_snapshot = 1 : no_snap_process          % Track positions
            
            c  = zeros( n_links , n_paths );          % The pattern coefficient matrix
            cp = zeros( n_links , n_paths );          % The phase coefficient matrix
            
            % Update the drifting angles, phases and delays.
            if drifting_precision > 0
                [ aoa, eoa, phase, delay(i_snapshot,:,:,:),...
                    aod_los, eod_los, aoa_los, eoa_los  ] =...
                    h_cb.update_drifting( i_snapshot );
            end
            
            % Include the direction on travel in the angles of arrival
            [ aoa_c, eoa_c, deg ] =...
                calc_rx_rotation( aoa, eoa, hdir(i_snapshot), gdir(i_snapshot) );
            
            % Apply the additional polarization rotation for the NLOS paths
            deg = deg - gamma( :,:,:, ones(1,size(aoa,4)) );

            switch drifting_precision
                case {0,1}      % Tx planar, Rx planar
                    
                    % Interpolate the receive antenna patterns
                    [Vr,Hr,Pr] = h_cb.par.rx_array(i_mobile).interpolate( aoa_c ,...
                        eoa_c, 1:n_rx, rx_azimuth_grid, rx_elevation_grid ,...
                        rx_patV, rx_patH, rx_element_pos );
                    Pr  = reshape( Pr,1,n_paths,n_rx );
                    
                    % Update the LOS Tx pattern
                    if drifting_precision == 1
                        [tmpVt, tmpHt, tmpPt] =...
                            h_cb.par.tx_array.interpolate( ...
                            aod_los, eod_los, 1:n_tx, tx_azimuth_grid, ...
                            tx_elevation_grid, tx_patV, tx_patH, tx_element_pos );
                        
                        Vt( 1,1,:,: ) = tmpVt( 1,o_subpaths,: );
                            Ht( 1,1,:,: ) = tmpHt( 1,o_subpaths,: );
                            Pt(1,1:n_clusters:n_paths,:)  = tmpPt( 1,o_subpaths,: );
                        if use_ground_reflection
                            Vt( 1,2,:,: ) = tmpVt( 1,o_subpaths*2,: );
                            Ht( 1,2,:,: ) = tmpHt( 1,o_subpaths*2,: );
                            Pt(1,2:n_clusters:n_paths,:)  = tmpPt( 1,o_subpaths*2,: );
                        end
                    end
                    
                otherwise % case {2,3,4,5}
                    
                    % Interpolate the receive antenna patterns (NLOS)
                    Vr  = zeros(1,n_clusters,n_subpaths,n_rx);
                    Hr  = zeros(1,n_clusters,n_subpaths,n_rx);
                    for i_rx = 1:n_rx
                        [Vr(:,:,:,i_rx),Hr(:,:,:,i_rx)] =...
                            h_cb.par.rx_array(i_mobile).interpolate( ...
                            aoa_c(:,:,:,i_rx) , eoa_c(:,:,:,i_rx) , i_rx , ...
                            rx_azimuth_grid , rx_elevation_grid , rx_patV , ...
                            rx_patH, rx_element_pos);
                    end
                    
                    if drifting_precision >= 4 % case {4,5}
                        
                        % Include the direction of travel in the angles of arrival
                        [ aoa_los_c, eoa_los_c, deg_LOS ] =...
                            calc_rx_rotation( aoa_los, eoa_los, hdir(i_snapshot), gdir(i_snapshot) );
                        
                        % Calculate the Rx-Array-Pattern LOS response for
                        % each element separately.
                        if use_ground_reflection
                            Vr_LOS  = zeros(1,2,1,n_rx,n_tx);
                            Hr_LOS  = zeros(1,2,1,n_rx,n_tx);
                            Vt_LOS  = zeros(1,2,1,n_rx,n_tx);
                            Ht_LOS  = zeros(1,2,1,n_rx,n_tx);
                        else
                            Vr_LOS  = zeros(1,1,1,n_rx,n_tx);
                            Hr_LOS  = zeros(1,1,1,n_rx,n_tx);
                            Vt_LOS  = zeros(1,1,1,n_rx,n_tx);
                            Ht_LOS  = zeros(1,1,1,n_rx,n_tx);
                        end
                        
                        for i_rx = 1:n_rx
                            [ Vr_LOS(1,:,1,i_rx,:), Hr_LOS(1,:,1,i_rx,:) ] =...
                                h_cb.par.rx_array(i_mobile).interpolate( ...
                                aoa_los_c(1,:,1,i_rx,:), eoa_los_c(1,:,1,i_rx,:), i_rx, ...
                                rx_azimuth_grid , rx_elevation_grid , rx_patV , ...
                                rx_patH, rx_element_pos);
                        end
                        
                        for i_tx = 1:n_tx
                            [ Vt_LOS(1,:,1,:,i_tx), Ht_LOS(1,:,1,:,i_tx) ] =...
                                h_cb.par.tx_array.interpolate( ...
                                aod_los(1,:,1,:,i_tx), eod_los(1,:,1,:,i_tx), i_tx, ...
                                tx_azimuth_grid, tx_elevation_grid, tx_patV,...
                                tx_patH, tx_element_pos );
                        end
                        
                        phase = reshape( phase, 1, n_paths, n_rx, n_tx );
                        
                    else % case {2,3}
                        
                        if use_ground_reflection
                            Vt_LOS  = zeros(1,2,1,n_rx,n_tx);
                            Ht_LOS  = zeros(1,2,1,n_rx,n_tx);
                            Pt_LOS  = zeros(1,2,1,n_rx,n_tx);
                        else
                            Vt_LOS  = zeros(1,1,1,n_rx,n_tx);
                            Ht_LOS  = zeros(1,1,1,n_rx,n_tx);
                            Pt_LOS  = zeros(1,1,1,n_rx,n_tx);
                        end
                        
                        % Calculate the LOS Tx pattern
                        [Vt_LOS(1,:,1,:,:), Ht_LOS(1,:,1,:,:), Pt_LOS(1,:,1,:,:)] =...
                            h_cb.par.tx_array.interpolate( ...
                            aod_los, eod_los, 1:n_tx, tx_azimuth_grid, ...
                            tx_elevation_grid, tx_patV, tx_patH, tx_element_pos );
                        
                        phase = reshape( phase, 1, n_paths, n_rx );

                    end
            end
            
            if drifting_precision == 0
                % Calculate the Doppler profile.
                doppler = reshape( cos(aoa_c+pi).*cos(eoa) ,1,n_paths );
            end

            switch use_polarization_rotation
                case 0 % Original WINNER model
                    xprmat = xprmat_winner;
                    
                case {1,2} % Only for drifting_precision = {0,1}
                    if drifting_precision <= 1
                        % Calculate a common XPRmat for all Tx and Rx antennas
                        xprmat(1,:) = cos(deg(:));
                        xprmat(2,:) = -sin(deg(:));
                        xprmat(3,:) = xprmat(2,:);
                        xprmat(4,:) = -xprmat(1,:);
                        
                        % Include circular phase offset
                        if use_polarization_rotation == 2
                            xprmat([1,2],:) = xprmat([1,2],:) .* conj(kappa([1,1],:));
                            xprmat([3,4],:) = xprmat([3,4],:) .* kappa([1,1],:);
                        end
                    end
                    
                case 3 % WINNER with LOS model
                    xprmat = diag([1,-1,1,-1]) * xprmat_winner;
            end
                        
            % The main loop to calculate the channel coefficients
            for i_rx = 1 : n_rx                               % Rx elements
                
                if drifting_precision >= 2
                    % Rx uses spherical waves. We need to calculate one
                    % XPRMAT for each Rx antenna separately.
                    switch use_polarization_rotation
                        case {1,2}
                            xprmat(1,:) = cos( reshape( deg(1,:,:,i_rx) ,1,n_paths ) );
                            xprmat(2,:) = -sin( reshape( deg(1,:,:,i_rx) ,1,n_paths ) );
                            xprmat(3,:) = xprmat(2,:);
                            xprmat(4,:) = -xprmat(1,:);
                            
                            % Include circular phase offset
                            if use_polarization_rotation == 2
                                xprmat([1,2],:) = xprmat([1,2],:) .* conj(kappa([1,1],:));
                                xprmat([3,4],:) = xprmat([3,4],:) .* kappa([1,1],:);
                            end
                    end
                end
                
                for i_tx = 1 : n_tx                           % Transmit elements
                    ind = (i_tx-1)*n_rx + i_rx;               % Index of element in c
                    
                    % Get the antenna pattern response for the NLOS paths
                    PatTx = [ reshape( Vt(1,:,:,i_tx) , 1,n_paths ) ;...
                        reshape( Ht(1,:,:,i_tx) , 1,n_paths ) ];
                    PatRx = [ reshape( Vr(1,:,:,i_rx) , 1,n_paths ) ;...
                        reshape( Hr(1,:,:,i_rx) , 1,n_paths ) ];
                    
                    % Update the LOS Tx pattern
                    if drifting_precision >= 2
                        los_ind = 1:n_clusters:n_paths;
                        PatTx(1,los_ind)    = Vt_LOS(1,1,1,i_rx,i_tx);
                        PatTx(2,los_ind)    = Ht_LOS(1,1,1,i_rx,i_tx);
                        if use_ground_reflection
                            gr_ind = 2:n_clusters:n_paths;
                            PatTx(1,gr_ind)    = Vt_LOS(1,2,1,i_rx,i_tx);
                            PatTx(2,gr_ind)    = Ht_LOS(1,2,1,i_rx,i_tx);
                        end
                    end
                    
                    % Update the LOS Rx pattern
                    if drifting_precision >= 4
                        PatRx(1,los_ind)    = Vr_LOS(1,1,1,i_rx,i_tx);
                        PatRx(2,los_ind)    = Hr_LOS(1,1,1,i_rx,i_tx);
                        if use_ground_reflection
                            PatRx(1,gr_ind)    = Vr_LOS(1,2,1,i_rx,i_tx);
                            PatRx(2,gr_ind)    = Hr_LOS(1,2,1,i_rx,i_tx);
                        end
                    elseif drifting_precision >= 2
                        Pt(1,los_ind,i_tx)  = Pt_LOS(1,1,1,i_rx,i_tx);
                        if use_ground_reflection
                            Pt(1,gr_ind,i_tx)  = Pt_LOS(1,2,1,i_rx,i_tx);
                        end
                    end
                    
                    if drifting_precision >= 4
                        % Update the LOS components of the XPRMAT for
                        % each Tx antenna. No update is required for mode
                        % {0,3}
                        switch use_polarization_rotation
                            case {1,2}
                                xprmat( 1,1:n_clusters:n_paths ) = cos( deg_LOS(1,1,1,i_rx,i_tx) );
                                xprmat( 2,1:n_clusters:n_paths ) = -sin( deg_LOS(1,1,1,i_rx,i_tx) );
                                xprmat( 3,1:n_clusters:n_paths ) = xprmat(2,1:n_clusters:n_paths);
                                xprmat( 4,1:n_clusters:n_paths ) = -xprmat(1,1:n_clusters:n_paths);
                        end
                    end

                    % Get the channel coefficients without random phases
                    c(ind,:) = sum( [ sum( PatTx .* xprmat([1 3],:)) ;...
                        sum( PatTx .* xprmat([2 4],:))] .* PatRx );
                    
                    % The phases
                    switch drifting_precision
                        case {0,1}
                            % In drifting mode, we have to update the coefficient
                            % matrix with the time-variant Doppler profile.
                            cp(ind,:) = exp( -1j*( pin +...
                                wave_no*( Pt(1,:,i_tx) + Pr(1,:,i_rx) ) +...
                                phase( 1 , : ) ));
                            
                        case {2,3}
                            % Tx phases come from the projection
                            % Rx phases are based on the values in "phase"
                            cp(ind,:) = exp( -1j*( pin +...
                                wave_no * Pt(1,:,i_tx) +...
                                phase( 1 , : , i_rx  )));
                            
                        case {4,5}
                            % The phases already contain the effect of the
                            % AoD. Hence, the parallel projection of the
                            % arrays is not needed.
                            cp(ind,:) = exp( -1j*( pin +...
                                phase( 1 , : , i_rx , i_tx )));

                    end
                end
            end
            
            % There atr be random phases in the sub-paths of the antenna
            % patterns. This changes the power when summing up the coefficients.
            
            % Combine antenna patterns and phases
            ccp = reshape( c.*cp, n_links , n_clusters , n_subpaths );
            ppat(:,:,i_snapshot) = sum( abs(ccp).^2, 3);        % Sum of the subpath powers
            cn(:,:,i_snapshot)   = sum( ccp, 3 );               % Complex sum
        end
        
        if drifting_precision == 0
            % Only one snapshot is calculated, the others are
            % emulated by phase rotation.
            
            % Combine pattern and phase for the first snapshopt
            c = c.*cp;
            
            for i_snapshot = 2 : n_snapshots
                
                % Generate rotating Dopplers for the sucessive snapshots
                cp = exp( -1j * wave_no * doppler * dist(i_snapshot) );
                cp = cp( ones(1,n_links) , : );
                
                % Combine antenna patterns and phases
                ccp = reshape( c.*cp, n_links , n_clusters , n_subpaths );
                ppat(:,:,i_snapshot) = sum( abs(ccp).^2, 3);    % Sum of the subpath powers
                cn(:,:,i_snapshot)   = sum( ccp, 3 );           % Complex sum
            end
        end
        
        % The path powers
        p_cl = h_cb.pow(i_mobile*ones(1,n_links),: );
        
        % The powers in the current channel coefficients (complex sum)
        p_coeff = mean( abs(cn).^2, 3 );
        
        % The powers of the antenna patterns at the given angles (power-sum)
        p_pat = mean( ppat,3 );
        
        % Correct the powers
        p_correct = sqrt( p_cl .* p_pat ./ p_coeff ./ n_subpaths );
        p_correct( p_pat < 1e-30 ) = 0; % Fix NaN caused by 0/0
        cn = p_correct(:,:,ones(1,n_snapshots)) .* cn;
        
        % Now we apply the K-Factor and the shadowing profile
        if drifting_precision > 0 && map_valid == 1
            
            % Get shadowing profile along the track from the correlation
            % map. The first vector is the K-Factor and the second vector
            % is the SF. The initial K-Factor is already applied in the
            % path powers. We thus need to correct for that factor.
            [sf,kf] = h_cb.par.get_sf_profile( ...
                h_cb.par.rx_track(i_mobile) , i_mobile );
            
            % Scaling factor for the KF
            kf  = kf./h_cb.par.kf(i_mobile);
            
            p1  = h_cb.pow(1);
            kf_power_scale = sqrt( 1+p1*(kf-1) );
            
            kf = sqrt(kf);
            
            % The path loss might be given manually together with the SF in
            % the  track object. In this case, we do not calculate it here
            if isempty( h_cb.par.rx_track(i_mobile).par ) || ...
                    isempty( h_cb.par.rx_track(i_mobile).par.pg )
                
                % Get the path loss
                [ path_loss , scale_sf ] = h_cb.par.get_pl( ...
                    h_cb.par.rx_track(i_mobile) , i_mobile );
            else
                % No path loss model is used when PL/SF are defined
                % manually.
                path_loss = 0;
                scale_sf = 1;
            end
            
            % We have the option to calculate the SF, PL and KF
            % antenna-dependent. This is activated when
            % simpar.drifting_precisioin is set to 3. However, SF
            % depends only on the rx position. We compute the values for
            % the tx-antennas here.
            if drifting_precision == 3
                rx_power = -path_loss + 10*log10( sf(:,:,o_tx) ) .* scale_sf;
            else
                rx_power = -path_loss + 10*log10( sf ) .* scale_sf;
            end
            rx_power = sqrt( 10.^( 0.1 * rx_power ) );
            
            if drifting_precision == 3
                kf = permute( reshape( kf(:,:,o_tx) , n_snapshots , [] ) , [2,3,1] );
                rx_power = rx_power ./ kf_power_scale(:,:,o_tx);
                rx_power = permute( reshape( rx_power , n_snapshots , [] ) , [2,3,1] );
            else
                o_tmp = ones(1,n_tx*n_rx);
                kf = permute( kf(:,o_tmp) , [2,3,1] );
                rx_power = rx_power ./ kf_power_scale;
                rx_power = permute( rx_power(:,o_tmp) , [2,3,1] );
            end
            
            cn(:,1,:) = cn(:,1,:).*kf;
            cn = cn.*rx_power(:,o_clusters,:);
            
        else
            % The path loss might be given manually together with the SF in
            % the  track object. In this case, we do not calculate it here
            if isempty( h_cb.par.rx_track(i_mobile).par ) || ...
                    isempty( h_cb.par.rx_track(i_mobile).par.pg )
                
                % The initial KF is already applied in path powers. Here,
                % we only need to apply the SF and the path loss.
                cn = cn * rx_power(i_mobile);
                
            else
                % Extract the path loss from the track object
                rx_power = h_cb.par.get_sf_profile( ...
                    h_cb.par.rx_track(i_mobile) , i_mobile );
                
                rx_power = sqrt( rx_power );
                
                o_tmp = ones(1,n_tx*n_rx);
                rx_power = permute( rx_power(:,o_tmp) , [2,3,4,1] );
                
                cn = cn.*rx_power(:,o_clusters,:);
            end
        end
        
        % Apply antenna coupling
        Ct = h_cb.par.tx_array.coupling;
        Cr = h_cb.par.rx_array(i_mobile).coupling.';
        
        % Reshape objects
        cn = reshape( cn , n_rx , n_tx , n_clusters , n_snapshots );
        
        % Apply the antenna coupling
        c = zeros( size(Cr,1) , size(Ct,2) , n_clusters , n_snapshots  );
        for i_snapshot = 1:n_snapshots
            for i_cluster = 1:n_clusters
                c(:,:,i_cluster,i_snapshot) = Cr * cn(:,:,i_cluster,i_snapshot) * Ct;
            end
        end
        clear cn
        
        if drifting_precision >= 2
            % When we use high precision, the delays on all elements are
            % different. However, antenna coupling will merge the
            % coefficients of different antennas. This needs to be
            % considered by the delays too.
            
            % The delays on different elements are weighted by the
            % powers in the coupling matrix.
            
            Cr_dl = zeros( size( Cr ));
            for i_rx = 1:size( Cr,1 )
                tmp = abs( Cr( i_rx , : ) ).^2;
                Cr_dl( i_rx , : ) = tmp./sum(tmp);
            end
            
            Ct_dl = zeros( size( Ct ));
            for i_tx = 1:size( Ct,2 )
                tmp = abs( Ct( : , i_tx ) );
                Ct_dl( : , i_tx ) = tmp./sum(tmp);
            end
            
            % Here, we scale the delays for each path by the
            % coupling.powers.
            
            delay = permute( delay, [3,4,2,1] );
            if drifting_precision < 4
                delay = delay( :,o_tx , : ,: );
            end
            dl = zeros( size(Cr,1) , size(Ct,2) , n_clusters , n_snapshots  );
            for i_snapshot = 1:n_snapshots
                for i_cluster = 1:n_clusters
                    dl(:,:,i_cluster,i_snapshot) = Cr_dl * delay(:,:,i_cluster,i_snapshot) * Ct_dl;
                end
            end
            h_channel(i_mobile) = channel( c , dl , initial_pos , false  );
        else
            h_channel(i_mobile) = channel( c , delay' , initial_pos , false  );
        end
        clear c
        
        h_channel(i_mobile).name = [ h_cb.name ,'_', h_cb.par.rx_track(i_mobile).name ];
        h_channel(i_mobile).rx_position = h_cb.par.rx_track(i_mobile).positions_abs;
        h_channel(i_mobile).tx_position = h_cb.par.tx_position;
        
        % Save Additional LSF and SSF information
        h_channel(i_mobile).par.ds_parset = h_cb.par.ds( i_mobile ); % [s]
        h_channel(i_mobile).par.kf_parset = 10*log10( h_cb.par.kf( i_mobile ) ); % [db]
        h_channel(i_mobile).par.pg_parset = 10*log10( mean(rx_power(:)).^2 ); % [db]
        h_channel(i_mobile).par.asD_parset = h_cb.par.asD( i_mobile ); % [deg]
        h_channel(i_mobile).par.asA_parset = h_cb.par.asA( i_mobile ); % [deg]
        h_channel(i_mobile).par.esD_parset = h_cb.par.esD( i_mobile ); % [deg]
        h_channel(i_mobile).par.esA_parset = h_cb.par.esA( i_mobile ); % [deg]
        
        h_channel(i_mobile).par.asD_cb = ...
            qf.calc_angular_spreads( h_cb.AoD( i_mobile,: ), h_cb.pow( i_mobile, : ) ) * 180/pi;
        h_channel(i_mobile).par.asA_cb = ...
            qf.calc_angular_spreads( h_cb.AoA( i_mobile,: ), h_cb.pow( i_mobile, : ) ) * 180/pi;
        h_channel(i_mobile).par.esD_cb = ...
            qf.calc_angular_spreads( h_cb.EoD( i_mobile,: ), h_cb.pow( i_mobile, : ) ) * 180/pi;
        h_channel(i_mobile).par.esA_cb = ...
            qf.calc_angular_spreads( h_cb.EoA( i_mobile,: ), h_cb.pow( i_mobile, : ) ) * 180/pi;
        
        h_channel(i_mobile).par.xpr_parset = 10*log10( h_cb.par.xpr( i_mobile ) ); % [db]
        h_channel(i_mobile).par.AoD_cb = h_cb.AoD( i_mobile,: ) * 180/pi; % [deg]
        h_channel(i_mobile).par.AoA_cb = h_cb.AoA( i_mobile,: ) * 180/pi; % [deg]
        h_channel(i_mobile).par.EoD_cb = h_cb.EoD( i_mobile,: ) * 180/pi; % [deg]
        h_channel(i_mobile).par.EoA_cb = h_cb.EoA( i_mobile,: ) * 180/pi; % [deg]
    end
    
end

if verbose && nargin == 1
    fprintf('] %5.0f seconds\n',round( etime(clock, tStart) ));
end

end
