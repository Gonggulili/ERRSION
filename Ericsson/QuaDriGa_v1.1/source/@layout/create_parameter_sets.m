function [ h_parset, h_cb ] = create_parameter_sets( h_layout, initialize, check_parfiles )
%CREATE_PARAMETER_SETS Creates parameter_set objects based on layout specification
%
%   [ h_parset , h_cb ] = CREATE_PARAMETER_SETS creates the
%   'parameter_set' objects and (optionally) the 'channel_builder' object.
%   This includes the parameter_set objects that handle correlated large
%   scale parameters, tracks and antenna arrays for the Rx and Tx side. The
%   parameter_set objects will be fully initialized and initial parameters
%   for the channel builder will be calculated.
%
%   Options:
%       initialize = 0 (default: 1)
%       Does not create parameter maps. This is useful when you want to
%       edit the parameters after creating the parameter_Set object.
%
%       check_parfiles = 0 (default: 1)
%       Disables the parsing of shortnames and the validity-check for the
%       config-files. This is useful, if you know that the parameters in
%       the files are valid. In this case, this saves some execution time.
%
%
% QuaDRiGa Copyright (C) 2011-2013 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
%
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% Parse Input variables
if exist( 'initialize' , 'var' ) && ~isempty( initialize )
    if ~( all(size(initialize) == [1 1]) ...
            && (isnumeric(initialize) || islogical(initialize)))
        error('??? "initialize" is invalid.')
    end
else
    initialize = 1;
end

% Parse Input variables
if exist( 'check_parfiles' , 'var' )
    if ~( all(size(check_parfiles) == [1 1]) ...
            && (isnumeric(check_parfiles) || islogical(check_parfiles)) ...
            && any( check_parfiles == [0 1] ) )
        error('??? "check_parfiles" must be 0 or 1')
    end
else
    check_parfiles = true;
end

scenarios   = h_layout.track(1).scenario(1,1);         	% Scenario of the first segment
tracks      = cell(1,h_layout.no_tx);
rx          = cell(1,h_layout.no_tx);
tx          = cell(1,h_layout.no_tx);
position    = cell(1,h_layout.no_tx);

% The following loop parses all tracks and extracts the scenarios

count = 1;
scenario_active = zeros(1,h_layout.no_tx);
for n = 1:h_layout.no_rx                                	% Do for each Rx
    
    % If we have less scenarios defined than transmitters in the layout
    if size( h_layout.track(n).scenario,1 ) == 1 && h_layout.no_tx > 1
        rx_scen = repmat( h_layout.track(n).scenario , h_layout.no_tx, 1 );
    else
        rx_scen = h_layout.track(n).scenario;
    end
    
    % We only process links where there is an enrty in the pairing matrix
    tx_selection = sort( h_layout.pairing(1, h_layout.pairing(2,:) == n ) );
    
    if isempty( tx_selection )
        disp(['Receiver "',h_layout.rx_name{n},'" has no transmitter.'])
    end
    
    if h_layout.track(n).no_segments == 1    % It the current track has only one segment
        subtrack = h_layout.track(n).copy;   % ... use the track.
        rename_track = false;
        
        % Check, if the initial position is on the track.
        % If it is not, set the initial position to the first snapshot of
        % the track. This is done automatically in "track.get_subtrack".
        if ~any( sum( subtrack.positions == 0 ) == 3 )
            sp = subtrack.positions( :, 1 );
            subtrack.initial_position = subtrack.initial_position + sp;
            for m=1:3
                subtrack.positions(m,:) = subtrack.positions(m,:) - sp(m);
            end
        end
    else                                % ... split the track in subtracks.
        subtrack = h_layout.track(n).get_subtrack;
        rename_track = true;
    end
    
    for m = 1:h_layout.track(n).no_segments              % Do for each segment
        if rename_track                             % Set name for subtrack
            subtrack(m).name = [h_layout.rx_name{n},'_seg',num2str(m,'%04u')];
        end
        
        pos = zeros(subtrack(m).no_snapshots,3);  	% Get all positions from current subtrack
        pos(:,1) = subtrack(m).positions(1,:) + subtrack(m).initial_position(1);
        pos(:,2) = subtrack(m).positions(2,:) + subtrack(m).initial_position(2);
        pos(:,3) = subtrack(m).positions(3,:) + subtrack(m).initial_position(3);
        
        % Each Rx can belong to several scenarios. One for each Tx.
        % o iterates over all tx
        for o = tx_selection
            % Check, if the scenario is already listed
            [ already_there , loc] = ismember( rx_scen(o,m) , scenarios );
            
            if ~already_there
                % If scenario does not exist, create it
                count = count + 1;
                scenarios(count) = rx_scen(o,m);
                position{ count,o }(:,1) = min( pos,[],1 )';
                position{ count,o }(:,2) = max( pos,[],1 )';
                scenario_active( count , o ) = 1;
                tracks{ count,o } = subtrack(m);
                rx{ count,o } = h_layout.rx_array(n);
                tx{ count,o } = h_layout.tx_array( o );
            else
                % Otherwise update position of the existing scenario
                position{ loc,o }(:,1) = min( [ pos ; position{loc,o}' ],[],1 )';
                position{ loc,o }(:,2) = max( [ pos ; position{loc,o}' ],[],1 )';
                scenario_active( loc , o ) = 1;
                tracks{ loc,o } = [ tracks{ loc,o },subtrack(m) ];
                rx{ loc,o } = [ rx{ loc,o },h_layout.rx_array(n) ];
                tx{ loc,o } = h_layout.tx_array( o );
            end
        end
    end
end

% The field names for the given parameters
par_fieldnames = {'ds','kf','pg','asD','asA','esD','esA','xpr'};

% Replace the scenario shortnames with their long form
[ sup_scenarios , file_names ] = ...
    parameter_set.supported_scenarios( check_parfiles );


% Add file names to list of supported scenarios
for i_scen = 1:numel( scenarios )
    if numel(scenarios{i_scen}) > 5 && ...
            ~isempty( regexp( scenarios{i_scen}(end-4:end) , '.conf', 'once' ) ) && ...
            regexp( scenarios{i_scen}(end-4:end) , '.conf' ) == 1
        file = dir(scenarios{i_scen});
        if ~isempty(file)
            sup_scenarios{end+1} = scenarios{i_scen};
            file_names{end+1} = [scenarios{i_scen},'.conf'];
        end
    end
end


% Create list of output scenarios
for n = 1:numel(scenarios)
    
    % Replace the scenario shortnames with their long form
    ind = strcmp( scenarios{n} , sup_scenarios );
    scenarios{n} = file_names{ind}(1:end-5);
    
    % Check the parameter files
    if check_parfiles
       parameter_set( scenarios{n},[],true );
    end
    
    for o = 1:h_layout.no_tx
        
        % When parameters are given in the tracks, we need to sort out, which
        % parameters belong to which transmitter. Hence, we need to split the
        % parameter-structures for different Txs here.
        
        if scenario_active(n,o)
            % Process parameters only when there are more then 2 Txs.
            
            for m = 1 : numel( tracks{n,o} )
                % Check each position
                
                if ~isempty( tracks{n,o}(m).par )
                    % Do if parameters are given (otherwise do nothing):
                    
                    not_copied = true;
                    for p = 1:numel( par_fieldnames )
                        tmp = tracks{n,o}(m).par.( par_fieldnames{p} );
                        if size( tmp,1 ) > 1
                            % Copy the track (since subtracks are only
                            % handles referring to the same object).
                            if not_copied
                                tracks{n,o}(m) = tracks{n,o}(m).copy;
                                not_copied = false;
                            end
                            % Use the row that matches the current
                            % tx-number.
                            tracks{n,o}(m).par.( par_fieldnames{p} ) = tmp(o,:);
                        end
                    end
                end
            end
        end
        
        % Create parameter set
        if o == 1
            h_parset(n,o) = parameter_set(scenarios{n},[],false);
        else
            h_parset(n,o) = h_parset(n,1).copy;
        end
        
        % Scenario-names should not contain underscores ('_')
        parset_name = regexprep( scenarios{n} , '_' , '-' );
        h_parset(n,o).name = [parset_name,'_',h_layout.tx_name{o}];
        
        h_parset(n,o).simpar = h_layout.simpar;
        h_parset(n,o).tx_position = h_layout.tx_position(:,o);
        
        if scenario_active(n,o)
            % Insert borders for map creation
            h_parset(n,o).positions = position{n,o};
            h_parset(n,o).map_extent =...
                [ position{n,o}(1,1) - h_parset(n,o).map_extension,...
                position{n,o}(1,2) + h_parset(n,o).map_extension;...
                position{n,o}(2,1) - h_parset(n,o).map_extension,...
                position{n,o}(2,2) + h_parset(n,o).map_extension];
            h_parset(n,o).samples_per_meter = h_layout.simpar.map_resolution;
            
            h_parset(n,o).rx_track   = tracks{n,o};
            h_parset(n,o).tx_array   = tx{ n,o };
            h_parset(n,o).rx_array   = rx{ n,o };
            
            % Insert segments for parameter generation
            h_parset(n,o).positions = cat( 2 , tracks{n,o}.initial_position );
        else
            h_parset(n,o).no_positions = 0;
        end
    end
end

switch initialize
    case 0  % Generate empty maps
        % Initialize parameters with default values and declare data and maps
        % as invalid
        h_parset.update_parameters(2);
        
    case 1
        % Initialize parameters with default values and declare data and maps
        % as invalid
        h_parset.update_parameters(2);
        
        % Get parameters either from provided LSPs in track.par or from maps
        h_parset.update_parameters(0);
end

if nargout == 2
    for n = 1:numel( h_parset )
        h_cb(n) = channel_builder( h_parset(n) );
    end
end

end

