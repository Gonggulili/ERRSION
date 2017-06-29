function update_parameters( h_parset, force )
%UPDATE_PARAMETERS Generates the  LSP maps and updates the parameters for all terminals
%
%   UPDATE_PARAMETERS calculates correlated large scale parameters for each
%   user position. Those parameters are needed by the channel builder class
%   to calculate initial parameters for each track or segment which are
%   then evolved into time varying channels.
%
%   By default, "update_parameters" reads the values given in the track
%   objects of the layout. If there are no values given or if parts of the
%   values are missing, the correlation maps are generated to extract the
%   missing parameters.
%
% Input:
%   force = 0 (default)
%   Tries to read the parameters from 'layout.track.par'. If they
%   are not provided or it they are incomplete, they are completed with
%   values from the LSP maps. If the maps are invalid (e.g. because
%   they have not been generated yet), new maps are  created.
%
%   force = 1
%   Creates new maps and reads the LSP from those maps. Values from
%   'layout.track.par' are ignored. Note that the parameters 'pg' and 'kf'
%   will still be taken from 'layout.track.par' when generating channel
%   coefficients.
%
%   force = 2
%   Creates dummy data for the maps and the LSPs. Any existing maps will be
%   deleted. Data and maps will be declared as invalid and the next time
%   when 'update_parameters' is called, new parameters are generated.
%   Values in 'layout.track.par' will NOT be affected.
%
%
% QuaDRiGa Copyright (C) 2011-2014 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
%
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% Parse input variables.
if nargin > 1
    force = real(force(1));
else
    force = 0;
end

% Update the correlations maps if necessary
if force < 0
    % Do nothing
    % When "update_parameters" is called recursively for each
    % "parameter_set"-object, we do not need to check this again.
    
elseif force == 1
    
    % Create new maps.
    % Ignore values from track in any case and use values from map
    h_parset = generate_correlation_maps( h_parset, 1 );
    force = -1;
    
elseif force == 2
    
    % Create empty maps.
    % Ignore values from track and use values from map.
    % Declare data as invalid.
    h_parset = generate_correlation_maps( h_parset, 0 );
    force = -1;
    
else % force == 0
    % Check if ALL parameters are provided by the tracks
    % (in this case, no maps are needed).
    
    par_complete = true;
    n=1; m=1;
    
    while par_complete && n <= numel(h_parset)
        
        % The loop stops when the first missing or incomplete set of
        % parameters is found. In this case, maps are generated.
        % When ALL parameters are complete (and the loop exits after
        % checking all), then no maps will be generated.
        
        par = h_parset(n).rx_track(m).par;
        if h_parset(n).no_positions == 0
            % There are no segments
            n=n+1;
            
        elseif isempty( par )
            % There are segments, but no parameters are given
            par_complete = false;
            
        elseif isempty( par.ds ) || isempty( par.kf ) || isempty( par.pg )  || ...
                isempty( par.asD ) || isempty( par.asA ) || isempty( par.esD ) || ...
                isempty( par.esA )
            % There are tracks with parameters, but the parameters are not
            % complete.
            
            par_complete = false;
        elseif m >= h_parset(n).no_positions
            % All segments of the current scenario have been checked
            n=n+1; m=1;
            
        else
            % There are unchecked segments in the current scenario.
            % Move to the next segment.
            m=m+1;
            
        end
    end
    
    % Calculate the maps only if there are parameters missing
    if ~par_complete
        % Determine which maps must be updated
        create_maps = ~cat( 2,h_parset.Pmap_valid ) & cat(2,h_parset.no_positions)>0;
        if any( create_maps )
            generate_correlation_maps( h_parset( create_maps ), 1 );
        end
        force = -2;     % Mixed values use maps for missing parameters
    else
        force = -3;     % Do not use maps and use only values from track
    end
    
end

% Parse Tx-Names from parameter_set names
tx_name = cell.empty;
for i_parset = 1:size(h_parset,2)
    % Parse tx index
    tmp = regexp( h_parset(1,i_parset).name , '_' );
    tx_name_local = h_parset(1,i_parset).name( tmp+1:end );
    tx_name{ i_parset } = tx_name_local;
end
no_tx = numel( tx_name );

% Iterate over all parsets
for i_parset = 1 : numel( h_parset )
    
    if h_parset(i_parset).no_positions > 0
        
        % Values of "force"
        %   -1    Ignore values from track in any case (force maps)
        %   -2    Mixed values use maps for missing parameters
        %   -3    Do not use maps and use only values from track (force track)
        
        data_valid = false;
        ksi = zeros(8,h_parset(i_parset).no_positions);
        
        if force ~= -3  % Interpolate values from the maps (for -1 and -2)
            
            if h_parset(i_parset).map_valid && ~h_parset(i_parset).data_valid
                % Maps are valid AND data is invalid
                
                if any( min( h_parset(i_parset).positions(1:2,:),[],2 ) ...
                        < h_parset(i_parset).map_extent(:, 1)) || ...
                        any( max( h_parset(i_parset).positions(1:2,:),[],2 ) ...
                        > h_parset(i_parset).map_extent(:, 2))
                    error('Some positions are outside the map.');
                end
                
                % Interpolate the maps using a custom 2D spline interpolation
                ksi(1:7,:) = spline_2d( h_parset(i_parset).positions(1:2,:) ,...
                    h_parset(i_parset).parameter_maps ,...
                    h_parset(i_parset).map_x_coord ,...
                    h_parset(i_parset).map_y_coord );
                
                % Transform to linear values
                ksi( [1,4:7],: ) = 10.^( ksi( [1,4:7],: ) );
                ksi( [2,3],: )   = 10.^( 0.1 * ksi( [2,3],: ) );
                
                % XPR
                log10_f_GHz  = log10( h_parset(i_parset).simpar.center_frequency / 1e9 );
                
                mu     = h_parset(i_parset).scenpar.xpr_mu;
                gamma  = h_parset(i_parset).scenpar.xpr_gamma;
                sigma  = h_parset(i_parset).scenpar.xpr_sigma;
                delta  = h_parset(i_parset).scenpar.xpr_delta;
                
                mu = mu + gamma * log10_f_GHz;
                sigma = sigma + delta * log10_f_GHz;
                sigma( sigma<0 ) = 0;
                
                xpr = randn(1,h_parset(i_parset).no_positions) * sigma + mu;
                ksi(8,:) = 10.^( 0.1*xpr );
                
                data_valid = true;
                
            elseif ~h_parset(i_parset).map_valid && ~h_parset(i_parset).data_valid
                % Maps are invalid AND data is invalid
                
                % Carrier frequency in GHz
                f_GHz = h_parset(i_parset).simpar.center_frequency / 1e9;
                
                % Get the average LSPs
                mu = [ h_parset(i_parset).scenpar.DS_mu,...
                    h_parset(i_parset).scenpar.KF_mu,...
                    0,...
                    h_parset(i_parset).scenpar.AS_D_mu,...
                    h_parset(i_parset).scenpar.AS_A_mu,...
                    h_parset(i_parset).scenpar.ES_D_mu,...
                    h_parset(i_parset).scenpar.ES_A_mu,...
                    h_parset(i_parset).scenpar.xpr_mu];
                
                gamma = [ h_parset(i_parset).scenpar.DS_gamma,...
                    h_parset(i_parset).scenpar.KF_gamma,...
                    0,...
                    h_parset(i_parset).scenpar.AS_D_gamma,...
                    h_parset(i_parset).scenpar.AS_A_gamma,...
                    h_parset(i_parset).scenpar.ES_D_gamma,...
                    h_parset(i_parset).scenpar.ES_A_gamma,...
                    h_parset(i_parset).scenpar.xpr_gamma];
                
                if any( gamma ~= 0 )
                    mu = mu + gamma * log10( f_GHz );
                end
                
                o = ones( 1,h_parset(i_parset).no_positions );
                for i_ksi = 1:8
                    switch i_ksi
                        case {2,8} % KF, XPR
                            ksi(i_ksi,:) = o.*10.^(0.1*mu(i_ksi));
                        case {3} % SF
                            ksi(i_ksi,:) = o;
                        otherwise % Ds and ASs
                            ksi(i_ksi,:) = 10.^mu(i_ksi);
                    end
                end
                data_valid = false;
                
            else % Data is valid
                
                ksi(1,:) = h_parset(i_parset).ds;
                ksi(2,:) = h_parset(i_parset).kf;
                ksi(3,:) = h_parset(i_parset).sf;
                ksi(4,:) = h_parset(i_parset).asD;
                ksi(5,:) = h_parset(i_parset).asA;
                ksi(6,:) = h_parset(i_parset).esD;
                ksi(7,:) = h_parset(i_parset).esA;
                ksi(8,:) = h_parset(i_parset).xpr;
                
                data_valid = true;
                
            end
            
        else
            ksi = zeros( 8,h_parset(i_parset).no_positions );
        end
        
        % Read values from the tracks and overwrite data from the maps
        if force ~= -1  % for force = -2 and -3
            
            % Parse Tx-Number from parameter_set name
            tmp = regexp( h_parset(i_parset).name , '_' );
            tx_name_local = h_parset(i_parset).name( tmp+1:end );
            tx_ind = strcmp( tx_name , tx_name_local );
            
            par_fieldnames = {'ds','kf','pg','asD','asA','esD','esA','xpr'};
            
            data_complete = true;
            for n = 1 : numel( h_parset(i_parset).rx_track )
                % Temporary copy of the par struct for faster access
                par = h_parset(i_parset).rx_track(n).par;
                
                if ~isempty( par )
                    seg_ind = h_parset(i_parset).rx_track(n).segment_index(end);
                   
                    for p = 1:8
                        % Temporarily read values
                        tmp = par.( par_fieldnames{p} );
                        
                        if size(tmp,1) == no_tx
                            t_ind = tx_ind;
                        elseif size(tmp,1) == 1
                            t_ind = 1;
                        elseif isempty( tmp )
                            % OK
                            data_complete = false;
                        else
                            error('??? Invalid dimensions of "track.par"');
                        end
                        
                        % Copy the data to the parameter matrix
                        if ~isempty( tmp )
                            if p==2 || p==3
                                ksi( p,n ) = 10.^( 0.1 *tmp(t_ind,seg_ind) );
                            elseif p == 8
                                ksi( p,n ) = 10.^( 0.1 *tmp(t_ind,end) );
                            else
                                ksi( p,n ) = tmp(t_ind,end);
                            end
                        end
                    end
                else
                    data_complete = false;
                end
            end
            
            % The data is valid if we either have a complete dataset or if
            % the data from the maps was valid.
            data_valid = data_valid | data_complete;
        end
        
        % Copy interpolated values to the output variables
        h_parset(i_parset).ds  = ksi(1,:);
        h_parset(i_parset).kf  = ksi(2,:);
        h_parset(i_parset).sf  = ksi(3,:);
        h_parset(i_parset).asD = ksi(4,:);
        h_parset(i_parset).asA = ksi(5,:);
        h_parset(i_parset).esD = ksi(6,:);
        h_parset(i_parset).esA = ksi(7,:);
        h_parset(i_parset).xpr = ksi(8,:);
        
        h_parset(i_parset).data_valid = data_valid;
    end
end

end