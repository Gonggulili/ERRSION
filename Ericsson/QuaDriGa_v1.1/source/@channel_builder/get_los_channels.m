function h_channel = get_los_channels( h_parset, precision, return_coeff, tx_array_mask )
%GET_LOS_CHANNELS Generates channel coefficients for the LOS path only
%
%   h_channel = GET_LOS_CHANNELS( h_parset ) generates static coefficients for the
%   LOS path only. This includes the following properties:
%
%     - antenna patterns
%     - polarization rotation for the LOS path
%     - plane-wave approximation of the phase
%     - path loss
%     - shadow fading
%
%    No further features of QuaDRiGa are used (i.e. no drifting, spherical
%    waves, time evolution, multipath fading etc.). This function can thus be
%    used to acquire quick previews of the propagation conditions for a given
%    layout.
%
%   precision = 'single';
%   The additional input parameter 'precision' can be used to set the
%   numeric precision to 'single', thus reducing the memory requirements
%   for certain computations.
%
%   return_coeff = 'coeff';
%   If this is set, only the raw channel coefficients are returned, but no
%   QuaDRiGa channel object is created. This may help to reduce memory
%   requirements.
%   In this case, 'h_channel' has the dimensions: [ n_rx, n_tx, n_pos ]
%
%   return_coeff = 'raw';
%   Same as 'coeff' but without applying the distance-dependant phase and
%   the path loss. The rx-antenna is assumed to be dual-polarized with two
%   elements (i.e. the rx interpolation is omitted). This mode is used
%   QuaDRiGa-internally by [array.combine_pattern].
%
%
% QuaDRiGa Copyright (C) 2011-2015 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
%
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

single_precision = true;
if ~exist('precision','var') || ~strcmp( precision , 'single' )
    precision = 'double';
    single_precision = false;
end

if exist('return_coeff','var') && ~isempty( return_coeff )
    switch return_coeff
        case 'coeff'
            return_coeff = true;
            raw_coeff = false;
        case 'raw'
            return_coeff = true;
            raw_coeff = true;
    end
else
    return_coeff = false;
    raw_coeff = false;
end

if return_coeff && numel( h_parset ) ~= 1
    error('Raw channel coefficients can only be generted for scalar parameter set opjects.')
else
    h_channel = channel.empty;
end

for i_parset = 1:numel( h_parset )
    
    if numel( unique(h_parset(i_parset).rx_array) ) > 1
        warning('MATLAB:channel_builder',...
            ['There is more than one Rx antenna arrays in "h_parset".\n ',...
            '"get_los_channels" can only use one array. Results might be erroneous.']);
    end
    
    lambda  = h_parset(i_parset).simpar.wavelength;
    n_positions = h_parset(i_parset).no_positions;
    wave_no = 2*pi/lambda;
    
    % Extract the travel directions at the initial position
    gdir = zeros(1,n_positions);
    if numel( h_parset(i_parset).rx_track ) == n_positions
        for i_pos = 1:n_positions
            if isempty( h_parset(i_parset).rx_track(i_pos).ground_direction )
                h_parset(i_parset).rx_track(i_pos).compute_directions;
            end
            [~,ind] = min( sum( h_parset(i_parset).rx_track(i_pos).positions.^2 ) );
            gdir(i_pos) = h_parset(i_parset).rx_track(i_pos).ground_direction( ind );
        end
    end
    
    % Get the arrival and departure angles
    angles = h_parset(i_parset).get_angles*pi/180;
    if single_precision
        angles = single( angles );
    end
    
    % Interpolate the patterns
    if exist( 'tx_array_mask','var' )
        [ Vt,Ht,Pt ] = h_parset(i_parset).tx_array.interpolate( angles(1,:) , angles(3,:), tx_array_mask );
        Ct = h_parset(i_parset).tx_array.coupling( tx_array_mask,tx_array_mask );
        n_tx = numel(tx_array_mask);
    else
        [ Vt,Ht,Pt ] = h_parset(i_parset).tx_array.interpolate( angles(1,:) , angles(3,:) );
        Ct = h_parset(i_parset).tx_array.coupling;
        n_tx = h_parset(i_parset).tx_array.no_elements;
    end
    
    if raw_coeff
        n_rx = 2;
    else
        [ Vr,Hr,Pr ] = h_parset(i_parset).rx_array(1).interpolate( angles(2,:)-gdir , angles(4,:) );
        n_rx = h_parset(i_parset).rx_array(1).no_elements;
    end
    
    % Calculate the distance-dependent phases
    if ~raw_coeff
        r = (h_parset(i_parset).positions(1,:) - h_parset(i_parset).tx_position(1)).^2 + ...
            (h_parset(i_parset).positions(2,:) - h_parset(i_parset).tx_position(2)).^2 + ...
            (h_parset(i_parset).positions(3,:) - h_parset(i_parset).tx_position(3)).^2;
        phase = 2*pi/lambda * mod( sqrt(r), lambda);
        
        if single_precision
            phase = single( phase );
        end
    end
    
    if single_precision
        wave_no = single( wave_no );
    end
    
    % Calculate the channel coefficients including polarization
    c = zeros( n_rx*n_tx , n_positions, precision );
    for i_tx = 1 : n_tx
        
        % Tx Patterns
        PatTx = [ reshape( Vt(1,:,i_tx) , 1,n_positions ) ;...
            reshape( Ht(1,:,i_tx) , 1,n_positions ) ];
        
        if raw_coeff
            
            % First component
            ind = (i_tx-1)*2 + 1;
            c(ind,:) = PatTx(1,:) .* exp( -1j*(  wave_no*( Pt(1,:,i_tx)  )));
            
            % Second component
            ind = ind + 1;
            c(ind,:) = PatTx(2,:) .* exp( -1j*(  wave_no*( Pt(1,:,i_tx)  )));
            
        else
            for i_rx = 1 : n_rx
                ind = (i_tx-1)*n_rx + i_rx;
                
                % Rx Patterns
                PatRx = [ reshape( Vr(1,:,i_rx) , 1,n_positions ) ;...
                    reshape( Hr(1,:,i_rx) , 1,n_positions ) ];
                
                % Coefficients and antenna-dependent phase offset
                c(ind,:) = ( PatTx(1,:) .* PatRx(1,:) - PatTx(2,:) .* PatRx(2,:) ).* ...
                    exp( -1j*(  wave_no*( Pt(1,:,i_tx) + Pr(1,:,i_rx) ) + phase ));
                
            end
        end
    end
    
    % Apply path loss
    if ~raw_coeff
        [ path_loss , scale_sf ] = h_parset(i_parset).get_pl;
        rx_power = -path_loss.' + 10*log10( h_parset(i_parset).sf ) .* scale_sf;
        rx_power = sqrt( 10.^( 0.1 * rx_power ) );
        c = c.*rx_power( ones(1,n_tx*n_rx) , : );
    end
    
    % Apply antenna coupling
    c = reshape( c , n_rx , n_tx , n_positions );
    
    if single_precision
        Ct = single( Ct );
        Cr = single( h_parset(i_parset).rx_array(1).coupling.' );
    else
        Cr = h_parset(i_parset).rx_array(1).coupling.';
    end
    
    if all(size(Ct) == [ n_tx , n_tx ]) && ...
            all(all( abs( Ct - eye(n_tx)) < 1e-10 )) && ...
            all(size(Cr) == [ n_rx , n_rx ]) && ...
            all(all( abs( Cr - eye(n_rx)) < 1e-10 ))
        
        % Both coupling matrixes are identity matrices.
        coeff = c;
        
    elseif all(size(Ct) == [ n_tx , n_tx ]) && ...
            all(all( abs( Ct - diag(diag(Ct)) ) < 1e-10 )) && ...
            all(size(Cr) == [ n_rx , n_rx ]) && ...
            all(all( abs( Cr - eye(n_rx)) < 1e-10 ))
        
        % The tx has a diagonal matrix and the rx an identity matrix
        coeff = zeros( n_rx , n_tx , n_positions , precision  );
        for i_tx = 1 : n_tx
            coeff( :,i_tx,: ) = c(:,i_tx,:) .* Ct( i_tx,i_tx );
        end
        
        
    elseif all(size(Cr) == [ n_rx , n_rx ]) && ...
            all(all( abs(Cr - eye(n_rx)) < 1e-10 ))
        
        % Only the Rx is an identity matrix.
        coeff = zeros( n_rx , size(Ct,2) , n_positions , precision  );
        for i_tx_in = 1 : n_tx
            for i_tx_out = 1 : size(Ct,2)
                coeff( :,i_tx_out,: ) = coeff( :,i_tx_out,: ) +...
                    c(:,i_tx_in,:) .* Ct(i_tx_in,i_tx_out );
            end
        end
        
    elseif all(size(Ct) == [ n_tx , n_tx ]) && ...
            all(all( abs(Ct - eye(n_tx)) < 1e-10 ))
        
        % Only the Tx is an identity matrix.
        coeff = zeros( size(Cr,1) , n_tx , n_positions , precision );
        for n = 1:n_positions
            coeff(:,:,n) = Cr * c(:,:,n);
        end
        
    else
        % Both coupling matrixes are not identity matrices.
        coeff = zeros( size(Cr,1) , size(Ct,2) , n_positions , precision );
        for n = 1:n_positions
            coeff(:,:,n) = Cr * c(:,:,n) * Ct;
        end
        
    end
    
    if return_coeff
        h_channel = coeff;
    else
        % Create output channel object
        coeff = reshape(coeff,size(coeff,1), size(coeff,2),1,n_positions);
        h_channel(i_parset) = channel( coeff , zeros(n_positions,1,precision) , 1 , false  );
        h_channel(i_parset).name = h_parset(i_parset).name;
        h_channel(i_parset).tx_position = h_parset(i_parset).tx_position;
        h_channel(i_parset).rx_position = h_parset(i_parset).positions;
    end
end

if numel( h_parset ) > 1
    h_channel = reshape( h_channel , size(h_parset) );
end

end