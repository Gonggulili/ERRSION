function h_parset = generate_correlation_maps( h_parset, initialize, vb_dots )
%GENERATE_CORRELATION_MAPS Generates a new set of correlation maps
%
%   GENERATE_CORRELATION_MAPS manually generates the correlation maps.
%   Correlation maps are needed to handle the distance-dependent correlation
%   between mobile terminals. E.g. when two terminals are close to each other,
%   they will see similar channels.
%
% QuaDRiGa Copyright (C) 2011-2016 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
%
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% Parse input parameters
if exist( 'initialize' , 'var' ) && ~isempty( initialize )
    if ~( all(size(initialize) == [1 1]) ...
            && (isnumeric(initialize) || islogical(initialize)) ...
            && any( initialize == [0 1] ) )
        error('??? "initialize" must be 0 or 1')
    end
else
    initialize = true;
end

verbose = h_parset(1).simpar.show_progress_bars;
if nargin < 3
    if verbose
        vb_dots = qf.init_progress_dots( ones(1,numel(h_parset)) );
    else
        vb_dots = zeros(1,numel(h_parset));
    end
end

if nargin < 3 && verbose && initialize
    fprintf('Parameters   [');
    tStart = clock;
end

if numel(h_parset) > 1
    for n=1:numel(h_parset)
        generate_correlation_maps( h_parset(n), initialize, vb_dots(n) );
    end
    
elseif h_parset.no_positions > 0
    
    if isempty( h_parset.scenpar )
       error('Scenario not set.');
    end
    
    if ~h_parset.LSP_matrix_isOK
        error('LSP_matrix is not positive-definite.');
    end

    % Carrier frequency in GHz
    f_GHz = h_parset.simpar.center_frequency / 1e9;
    
    % Get the average LSPs
    mu = [ h_parset.scenpar.DS_mu,...
        h_parset.scenpar.KF_mu,...
        0,...
        h_parset.scenpar.AS_D_mu,...
        h_parset.scenpar.AS_A_mu,...
        h_parset.scenpar.ES_D_mu,...
        h_parset.scenpar.ES_A_mu];
    
    gamma = [ h_parset.scenpar.DS_gamma,...
        h_parset.scenpar.KF_gamma,...
        0,...
        h_parset.scenpar.AS_D_gamma,...
        h_parset.scenpar.AS_A_gamma,...
        h_parset.scenpar.ES_D_gamma,...
        h_parset.scenpar.ES_A_gamma];
    
    if any( gamma ~= 0 )
        mu = mu + gamma * log10( f_GHz );
    end

    % Get the std. of the LSPs
    sigma = [ h_parset.scenpar.DS_sigma,...
        h_parset.scenpar.KF_sigma,...
        h_parset.scenpar.SF_sigma,...
        h_parset.scenpar.AS_D_sigma,...
        h_parset.scenpar.AS_A_sigma,...
        h_parset.scenpar.ES_D_sigma,...
        h_parset.scenpar.ES_A_sigma];
    
    delta = [ h_parset.scenpar.DS_delta,...
        h_parset.scenpar.KF_delta,...
        h_parset.scenpar.SF_delta,...
        h_parset.scenpar.AS_D_delta,...
        h_parset.scenpar.AS_A_delta,...
        h_parset.scenpar.ES_D_delta,...
        h_parset.scenpar.ES_A_delta];
    
    if any( delta ~= 0 )
        sigma = sigma + delta * log10( f_GHz );
        sigma( sigma<0 ) = 0;
    end
    has_sigma = sigma > 0;
    
    % Get the autocorrelation distances and scale them to the map resolution.
    delta = [ h_parset.scenpar.DS_lambda,...
        h_parset.scenpar.KF_lambda,...
        h_parset.scenpar.SF_lambda,...
        h_parset.scenpar.AS_D_lambda,...
        h_parset.scenpar.AS_A_lambda,...
        h_parset.scenpar.ES_D_lambda,...
        h_parset.scenpar.ES_A_lambda];
    
    samples_per_meter = h_parset.samples_per_meter;
    samples_per_meter_reasonable = 4 / min( delta( delta>0 ) );
    
    if samples_per_meter == 0
        % Set the sampling rater automatically
        samples_per_meter = samples_per_meter_reasonable;
        
    elseif samples_per_meter < samples_per_meter_reasonable
        warning('QuaDRiGa:parameter_set:generate_correlation_maps',...
            ['Map resolution too low. Suggested value: ',...
            num2str( samples_per_meter_reasonable ),' samples / m ']);
    end
    
    h_parset.map_extension = ceil( 5 / samples_per_meter );
    
    delta = delta .* samples_per_meter;
        
    % Special case:
    % If the elevation of departure spread is distance-dependent, we need
    % to generate a map anyway, even if the ES_D_sigma is 0.
    if h_parset.scenpar.ES_D_mu_A ~= 0
        has_sigma( 6 ) = true;
    end
    
    % Extract the minimum and maximum positions from the list of positions
    min_pos = min( h_parset.positions-h_parset.map_extension,[],2 );
    max_pos = max( h_parset.positions+h_parset.map_extension,[],2 );
    
    % Determine the map edges taking an existing grid into account
    x_min = min( [ h_parset.map_extent(1, 1) , min_pos(1) ] );
    x_max = max( [ h_parset.map_extent(1, 2) , max_pos(1) ] );
    y_min = min( [ h_parset.map_extent(2, 1) , min_pos(2) ] );
    y_max = max( [ h_parset.map_extent(2, 2) , max_pos(2) ] );
    
    % Extending the map while keeping the sampling rate constant
    h_parset.map_extent = [x_min, x_max; y_min, y_max];
    
    % Set the samples per meter. This automatically adjusts the map size.
    h_parset.samples_per_meter = samples_per_meter;
    
    % Get map dimensions
    no_map_y = h_parset.map_size(2);
    no_map_x = h_parset.map_size(1);
    
    if initialize
        
        switch h_parset.simpar.use_map_algorithm
            case 1
                h_parset = generate_map_winner( h_parset, delta, has_sigma, vb_dots );
            case 2
                h_parset = generate_map_filter( h_parset, delta, has_sigma, vb_dots );
        end
        
        % Copy map data
        names = {'ds_map','kf_map','sf_map','asD_map','asA_map','esD_map','esA_map'};
        ind = find( has_sigma );
        maps = zeros( h_parset.map_size(2)*h_parset.map_size(1), sum(has_sigma) , 'single' );
        for n = 1 : sum(has_sigma)
            maps(:,n) = h_parset.( names{ind(n)} )(:);
        end
        
        % If there is no variance, no inter-parameter correlation needs to
        % be applied. 
        R_sqrt = sqrtm( h_parset.LSP_xcorr_matrix );
        R_sqrt = R_sqrt( has_sigma,has_sigma );
        
        % Apply inter-parameter correlation
        maps = maps * R_sqrt;
        
        % Reapply the variances after applying the correlation values
        mu_map  = mean( maps );
        sig_map = std( maps );
        for n = 1 : sum(has_sigma)
            maps(:,n) = (maps(:,n)-mu_map(n))./sig_map(n);
        end

        % Transform Normal distributed maps to scenario specific distributions
        % Copy correlated data back to the initial function
        cnt = 1;
        for n = 1:7
            a = mu(n);
            if n == 6 && h_parset.scenpar.ES_D_mu_A ~= 0
                a = 0;
            end
            if has_sigma(n)
                h_parset.(names{n}) =...
                    reshape( sigma(n)*maps(:,cnt) + a, no_map_y, no_map_x );
                cnt = cnt + 1;
            else
                h_parset.(names{n}) = a;
            end
        end
        
        % Apply the distance-dependent ESD
        if h_parset.scenpar.ES_D_mu_A ~= 0
            x = ones( h_parset.map_size(2) , 1) *...
                ( single( h_parset.map_x_coord - h_parset.tx_position(1) ) ./ 1000 );
            y = ( single( h_parset.map_y_coord - h_parset.tx_position(2)  ) ./ 1000 ).' * ...
                ones( 1, h_parset.map_size(1) ) ;
            d_2d_km = sqrt( x.^2 + y.^2 );
            clear x y
            
            esd_mu = h_parset.scenpar.ES_D_mu_A .* d_2d_km + mu(6);
            clear d_2d_km
            
            esd_mu( esd_mu < h_parset.scenpar.ES_D_mu_min ) = h_parset.scenpar.ES_D_mu_min;
            h_parset.esD_map = h_parset.esD_map + esd_mu;
            clear esd_mu
        end
        
        h_parset.Pmap_valid = true;
        
    else
        % Initialize the maps with default parameters
        h_parset.Pmap_ds  = mu(1);
        h_parset.Pmap_kf  = mu(2);
        h_parset.Pmap_sf  = 0;
        h_parset.Pmap_asD = mu(4);
        h_parset.Pmap_asA = mu(5);
        h_parset.Pmap_esD = mu(6);
        h_parset.Pmap_esA = mu(7);
        
        h_parset.Pmap_valid = false;
    end
    h_parset.Pdata_valid = false;
    
else
    if verbose && initialize
        for n = 1:vb_dots
            fprintf('o');
        end
    end
end

if nargin < 3 && verbose && initialize
    fprintf('] %5.0f seconds\n',round( etime(clock, tStart) ));
end

end
