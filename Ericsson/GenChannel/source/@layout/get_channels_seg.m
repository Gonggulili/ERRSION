function [ h_channel, h_cb ] = get_channels_seg( h_layout, tx, rx, seg, overlap )
%GET_CHANNELS_SEG Returns the channel coefficients for a single segment only
%
%   This function can be used to obtain the channel coefficients for a
%   single segment or single rx-tx combination only. Thus, the channel
%   model can be run in "streaming-mode", where updates are provided on the
%   fly. This can significantly reduce the memory requirements for long
%   time-sequences. 
%   (However, the parameter maps still need to be generated for the entire
%   scenario.)
%
%   Input:
%       tx      
%           The index of the transmitter (e.g. the BS)
%
%       rx      
%           The index of the receiver, or track (e.g. the BS)
%       
%       seg     
%           The segment indices on the the track. If it is not provided or
%           empty, the entire track is returned. It is also possible to
%           concatenate successive segments, i.e.: [1:3] or [3:5], etc.
%
%       overlap   
%           The overlapping fraction for the channel merger. 
%           Default is 0.5
%
%   Output:
%       h_channel   
%           The channel for the segment.
%
%       h_cb    
%           The channel builder object with initialized parameters
%           (Warning: Do not edit the parameters, this might lead to errors
%           or wrong results. Use customized config files instead.)
%
%   Features:
%       - Parameter maps will be deleted after the parameters were
%         extracted to free memory before creating the channels.
%       - Caching will be used to avoid multiple calculation of the same
%         overlapping regions.
%       - Preinitialization will be used to return the same coefficients
%         for successive calls.
%
% QuaDRiGa Copyright (C) 2011-2015 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
%
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

persistent h_cb_int chksum ind chan_cache chan_cache_ind

% Parse the input variables
if ~exist( 'tx' , 'var' ) || isempty( tx ) || numel(tx) > 1
    error('??? "tx" must be given as scalar integer number')
end

if ~exist( 'rx' , 'var' ) || isempty( rx ) || numel(rx) > 1
    error('??? "rx" must be given as scalar integer number')
end

if exist( 'seg' , 'var' ) && ~isempty( seg )
    seg = seg( : ).';
    if ~all( diff(seg) == 1 )
        error('??? "seg" must be a consecutive list of segments')
    end
else
    seg = [];
end

if exist( 'overlap' , 'var' ) && ~isempty( overlap )
    if ~( isnumeric(overlap) && all(size(overlap) == [1 1]) && isreal(overlap) ...
            && overlap<=1 && overlap>=0 )
        error('??? "overlap" must be scalar, and in between 0 and 1')
    end
else
    overlap = 0.5;
end

% Get a checksum indicating if the layout-object changed
X = h_layout.no_tx + h_layout.no_rx + overlap;
for ii = 1 : h_layout.no_tx
    X = X + sum( cast( h_layout.tx_name{ii} , 'double' ) );
end
for ii = 1 : h_layout.no_rx
    % Track names
    X = X + sum( cast( h_layout.track(ii).name , 'double' ) );
    
    % Initial Positions
    X = X + sum( h_layout.track(ii).initial_position(:) );
    
    % No. snapshots
    X = X + sum( h_layout.track(ii).no_snapshots );
    
    % Positions
    X = X + sum( h_layout.track(ii).positions(:) );
    
    % Segment index
    X = X + sum( h_layout.track(ii).segment_index(:) );
    
    % Scenarios
    for ij = 1 : numel( h_layout.track(ii).scenario )
        X = X + sum( cast( h_layout.track(ii).scenario{ij},'double' ) );
    end
end

if isempty( chksum ) || X ~= chksum
    chksum = X;
    update = true;
else
    update = false;
end

% Create correlated parameters. This only needs to be done once. If it has
% been done, we do not need to check it again
if isempty( h_cb_int ) || update
    
    % Disable paring mode
    h_layout.set_pairing;
    
    % Generate the parameters and discard the maps
    h_layout.generate_parameters( overlap, 1 );
    
    % Create the parameter sets. Those do not hold the maps, but only a
    % copy of the parameters from the layout.
    [ ~, h_cb_int ] = h_layout.create_parameter_sets;
    
    % Initialize the parameters for the channel builder.
    % That makes sure that successive calls of "get_channels_seg" return
    % the same output.
    for i_cb = 1:numel(h_cb_int)
        h_cb_int( i_cb ).init_parameters;
        
        % Set the LOS phase offset to 0.
        % This is needed since otherwise the merger would track the phase
        % offset.
        n_path = size(h_cb_int( i_cb ).taus,2);
        h_cb_int( i_cb ).pin(:,1:n_path:end) = 0;
    end
    
    % Create an index linking the contents of the layout to the parsets
    % The rows of the index are as follows:
    %   1   Tx Number
    %   2   Rx Number
    %   3   Segment Number
    %   4   Column number of the previous segment
    %   5   Column number of the following segment
    %   6   Index of the parameter_set object
    %   7   Position number in the parameter_set object
    
    index = [];
    for i_rx = 1:h_layout.no_rx
        
        seg_ind = 1 : h_layout.track( i_rx ).no_segments;
        rx_ind  = ones( 1,numel( seg_ind ) ) * i_rx;
        
        prev = (seg_ind - 1) + size( index , 2 );
        next = (seg_ind + 1) + size( index , 2 );
        if h_layout.track( i_rx ).closed == 1
            prev(1) = prev(end) + 1;
            next(end) = prev(2);
        else
            prev(1) = 0;
            next(end) = 0;
        end
        
        index = [index, [ rx_ind; seg_ind; prev; next ]];
    end
    
    % Set the first Tx
    index = [ ones(1,size(index,2));index ];
    
    % Set the remaining Tx
    tmp = index;
    for i_tx = 2 : h_layout.no_tx
        prev = tmp( 4,: ) + size( index , 2 );
        prev( tmp( 4,: )==0 ) = 0;
        
        next = tmp( 5,: ) + size( index , 2 );
        next( tmp( 5,: )==0 ) = 0;
        
        index = [index, [ tmp(1,:)*i_tx; tmp(2,:); tmp(3,:); prev; next ]];
    end
    
    % Calculate the mapping of par_set to tx
    partx = zeros( numel(h_cb_int), 1 );
    parrx = {};
    for i_par = 1:numel( h_cb_int )
        for i_tx = 1:h_layout.no_tx
            tmp = strfind( h_cb_int( i_par ).name , h_layout.tx_name{ i_tx } );
            if ~isempty( tmp )
                partx( i_par ) = i_tx;
            end
        end
        
        no_pos = h_cb_int( i_par ).par.no_positions;
        
        parrx{ i_par } = zeros( 2,no_pos  );
        for i_pos = 1 : no_pos
            tmp = h_cb_int( i_par ).par.rx_track( i_pos ).name;
            ct  = regexp( tmp,'_' );
            if isempty( ct )
                parrx{ i_par }(1,i_pos) = find(  strcmp( tmp , h_layout.rx_name ) );
                parrx{ i_par }(2,i_pos) = 1;
            else
                parrx{ i_par }(1,i_pos) = find( strcmp( tmp(1:ct-1) , h_layout.rx_name ) );
                parrx{ i_par }(2,i_pos) = str2double( tmp(ct+4:end) );
            end
        end
    end
    
    % Index the parameter_set
    par_ind = zeros( 1,size(index,2) );
    seg_ind = zeros( 1,size(index,2) );
    for ii = 1:size(index,2)
        % The list of possible parsets matching the tx
        pari = find( partx == index( 1,ii ) );
        
        for ij = 1:numel( pari )
            tmp = parrx{ pari(ij) };
            rxi = tmp(1,:) == index( 2,ii ) & tmp(2,:) == index( 3,ii );
            if any( rxi )
                par_ind(ii) = pari( ij );
                seg_ind(ii) = find( rxi );
            end
        end
    end
    
    ind = [ index ; par_ind ; seg_ind ];
end

% Pass channel builder object to the outside
h_cb = h_cb_int;

% Test if the tx and rx index are there
tx_rx_ind = ind(1,:) == tx & ind(2,:) == rx;
if ~any( tx_rx_ind )
    error('??? Transmitter or receiver is not in the layout.')
end

% If seg is empty, return all segments of the track
if isempty( seg )
    seg = ind( 3, tx_rx_ind );
    use_all_seg = true;
else
    use_all_seg = false;
end

% Temporary disable the progress bar
show_progress_bars = h_layout.simpar.show_progress_bars;
h_layout.simpar.show_progress_bars = 0;

chan_curr = channel.empty;
for i_seg = 1 : numel( seg )
    
    % Get the segment index
    ii = find( ind(1,:) == tx & ind(2,:) == rx & ind(3,:) == seg(i_seg) );
    if numel( ii ) ~= 1
        error(['??? Segment index ',num2str( seg(i_seg) ),' not found in track.'])
    end
    
    % Get the current channel segment
    if i_seg == 1 && ~isempty( chan_cache_ind ) && chan_cache_ind == ii
        chan_curr(1) = chan_cache;
    else
        chan_curr(i_seg) = get_channels_seg_unmerged( h_cb, ind(7,ii), ind(6,ii) );
    end
    
    if i_seg == numel(seg)
        
        if ind(5,ii) ~= 0
            % There is a next segment that needs to be merged
            
            % Get the channel of the following segment
            ij = ind(5,ii);
            chan_next = get_channels_seg_unmerged( h_cb, ind(7,ij), ind(6,ij) );
            
            % Special case: The next segment is the first one in case of
            % closed tracks. The segment number must be changed in order to
            % get the correct merger behavior.
            if ij == 1
                tmp = regexp( chan_next.name, '_seg' );
                chan_next.name = [ chan_next.name(1:tmp),'seg',...
                    num2str( ind(3,ii)+1 , '%04d' ) ];
            end
            
            % Merge the channels
            c = [ chan_curr , chan_next ];
            h_channel = c.merge(overlap);
            
            % Remove parts of the next segment that fall into the merger
            nv = 1 : h_channel.no_snap - chan_next.no_snap +...
                chan_next.initial_position - 1;
            
            delay = h_channel.delay;
            h_channel.coeff = h_channel.coeff(:,:,:,nv);
            if h_channel.individual_delays
                h_channel.delay = delay(:,:,:,nv);
            else
                h_channel.delay = delay(:,nv);
            end
            
            % Store the next channel object in cache in case we need it in the next
            % iteration.
            if ij ~= 1
                chan_cache      = chan_next;
                chan_cache_ind  = ij;
            else
                % Clear Cache
                chan_cache      = [];
                chan_cache_ind  = [];
            end
        else
            % There is no next segment
            h_channel = chan_curr.merge;

            % Clear Cache
            chan_cache      = [];
            chan_cache_ind  = [];
        end
        
        % Set the name
        if ~use_all_seg
            tmp = [h_channel.name , '_Seg',num2str(seg,'%04d+')];
            h_channel.name = tmp(1:end-1);
        end
        
    end
end

% Set progress-bar state to old value
h_layout.simpar.show_progress_bars = show_progress_bars;

end

function  chan = get_channels_seg_unmerged( h_cb, is, ip )

% Create temporary parameter set object for the current channel
h_par               = parameter_set( h_cb( ip ).par.scenario, [] , false );
h_par.name          = h_cb( ip ).par.name;
h_par.simpar        = h_cb( ip ).par.simpar;
h_par.tx_array      = h_cb( ip ).par.tx_array;
h_par.rx_array      = h_cb( ip ).par.rx_array( is );
h_par.rx_track      = h_cb( ip ).par.rx_track( is );
h_par.positions     = h_cb( ip ).par.positions(:,is);
h_par.tx_position   = h_cb( ip ).par.tx_position;
h_par.ds            = h_cb( ip ).par.ds( is );
h_par.kf            = h_cb( ip ).par.kf( is );
h_par.sf            = h_cb( ip ).par.sf( is );
h_par.asD           = h_cb( ip ).par.asD( is );
h_par.asA           = h_cb( ip ).par.asA( is );
h_par.esD           = h_cb( ip ).par.esD( is );
h_par.esA           = h_cb( ip ).par.esA( is );
h_par.xpr           = h_cb( ip ).par.xpr( is );
h_par.map_extent    = h_cb( ip ).par.map_extent;
h_par.data_valid    = true;

h_cb_tmp            = channel_builder( h_par );
h_cb_tmp.taus       = h_cb( ip ).taus( is,: );
h_cb_tmp.pow        = h_cb( ip ).pow( is,: );
h_cb_tmp.AoD        = h_cb( ip ).AoD( is,: );
h_cb_tmp.AoA        = h_cb( ip ).AoA( is,: );
h_cb_tmp.EoD        = h_cb( ip ).EoD( is,: );
h_cb_tmp.EoA        = h_cb( ip ).EoA( is,: );
h_cb_tmp.xpr        = h_cb( ip ).xpr( is,:,: );
h_cb_tmp.pin        = h_cb( ip ).pin( is,: );
h_cb_tmp.kappa      = h_cb( ip ).kappa( is,:,: );
h_cb_tmp.random_pol = h_cb( ip ).random_pol( :,:,is );
h_cb_tmp.subpath_coupling = h_cb( ip ).subpath_coupling(:,:,:,is);

chan = h_cb_tmp.get_channels;

end
