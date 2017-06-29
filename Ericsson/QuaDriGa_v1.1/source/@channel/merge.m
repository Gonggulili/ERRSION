function c = merge( h_channel, overlap, optimize, verbose )
%MERGE Combines channel segments into a continuous time evolution channel
%
%   c = MERGE( overlap ) Combines channel segments into a continuous time
%   evolution channel. The optional parameter 'overlap' scales the time
%   span during which the merging process takes place. It can have values
%   in between 0 and 1. A value of 0 disables the merging process and the
%   channel segments are simply concatenated. A value of 1 constantly
%   merges the channels. The default setting is 0.5.
%
%   c = MERGE( overlap, optimize ) The channel merger tries to
%   automatically optimize the pairing of the taps (i.e. one tap if the old
%   segment ramps down and one of the new ramps up). This is enabled by
%   default, but it is computing intensive. For quicker results, it can be
%   disabled by setting 'optimize' to 0.
%
%   c = MERGE( overlap, optimize , verbose ) Shows a progress bar
%   (default). This can be disabled by setting 'verbose' to 0.
%
% The channel merger implements the continuous time evolution with smooth
% transitions between segments. Each segment of a track is split in two
% parts: an overlapping area with the previous segment and an exclusive
% part with no overlapping. Each segment is generated independently by the
% channel builder. However, the distance dependent autocorrelation of the
% large scale parameters was considered when the parameters were drawn from
% the corresponding statistics.
%
% Transition from segment to segment is carried out by replacing taps of
% the previous segment by the taps of the current segment, one by one. The
% modeling of the birth/death process is done as published in the
% documentation of the WIM2 channel model. The route between adjacent
% channel segments is split into sub-intervals equal to the minimum number
% of taps in both overlapping segments. During each sub-interval the power
% of one old tap ramps down and one new tap ramps up. Power ramps are
% modeled by a modified sinus function to allow smooth transitions.
%
%   0.5*(1+sin((x-0.5)*pi))
%
% where x is the running index of the ramping window. x must be in between
% 0 and 1. Taps from the old and new segments are coupled based on their
% power. If the number of clusters is different in the channel segments,
% the weakest clusters are ramped up or down without a counterpart from the
% new/old segment. The merging is only done for the NLOS components since
% the LOS component has a deterministic behavior. The LOS component is thus
% just scaled in power.
%
% REFERENCE
% The time evolution is implemented as described in the WIM2
% documentation. "Kyösti, P.; Meinilä, J.; Hentilä, L. & others;
% {IST-4-027756 WINNER II D1.1.2 v.1.1}: WINNER II Channel Models; 2007".
%
% QuaDRiGa Copyright (C) 2011-2013 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
%
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% Parse input arguments
if exist( 'overlap' , 'var' )
    if ~( isnumeric(overlap) && all(size(overlap) == [1 1]) && isreal(overlap) ...
            && overlap<=1 && overlap>=0 )
        error('??? "overlap" must be scalar, and in between 0 and 1')
    end
else
    overlap = 0.5;
end

if exist( 'optimize' , 'var' )
    if ~( isnumeric(optimize) || ~isnan(optimize) || all(size(optimize) == [1 1]) )
        error('"optimize" must be numeric.')
    else
        optimize = logical(optimize);
    end
else
    optimize = true;
end

if exist( 'verbose' , 'var' )
    if ~( isnumeric(verbose) || ~isnan(verbose) || all(size(verbose) == [1 1]) )
        error('"verbose" must be numeric.')
    else
        verbose = logical(verbose);
    end
else
    verbose = true;
end

% Get the number of channels
n_channel = numel(h_channel);

% Show a progress bar if there are more than 10 segments
if n_channel < 10
    verbose = false;
end
if verbose
    fprintf('Merging      [');
    vb_dots = 50;
    tStart = clock;
    m0=0;
end

id_trk = cell(n_channel,1);                 % Reserve memory
id_seg = cell(n_channel,1);

% store the QuaDRiGa version number
qv = h_channel(1).version;

% Parse the name stings from the channel objects. This finds matching
% segments and combines them.

only_one_segment = true;
for i_channel = 1 : n_channel               % Do for each channel
    name = h_channel(i_channel).name;       % read channel name
    tmp = regexp( name , '_' );             % split "Scen_Tx_Rx_Seg"
    
    if numel( tmp ) == 1                    % We have already "Tx_Rx"
        tx = name(1:tmp(1)-1);              % store tx name
        rx = name(tmp(1)+1:end);            % store rx name
        seg = 'seg0001';                    % set segment number to 1
    else
        tx = name(tmp(1)+1:tmp(2)-1);       % store tx name
        if numel(tmp)==2                    % if there is only one segment ...
            rx = name(tmp(2)+1:end);        % store rx name
            seg = 'seg0001';                % set segment number to 1
        else
            rx = name(tmp(2)+1:tmp(3)-1);   % store rx name
            seg = name(tmp(3)+1:end);       % store segment name
            only_one_segment = false;
        end
    end
    
    id_trk{i_channel} = [tx,'_',rx];        % create a unique track id ...
    id_seg{i_channel} = [tx,'_',rx,'_',seg];   % ... and a segment id
end
[~,ind] = sort( id_seg );                   % sort names by segment

% Sort the channels
h_channel = h_channel( ind );
id_trk = id_trk( ind );

% Get the unique tracks
% The following code implements "unique(id_trk)" but does not sort the
% values in id_trk.

trk = id_trk;
if ~only_one_segment
    pt = 1;
    while pt < numel(trk)
        ii = find( strcmp( trk(pt) , trk(pt+1:end) ) )+pt;
        if ~isempty( ii )
            trk = trk( setdiff( 1:numel(trk), ii ) );
        end
        pt  = pt + 1;
    end
end
n_trk = numel(trk);

% The channel object can have an additional field for the position data.
% However, this field is mandatory. If it is empty, no processing of the
% positions will be done.

process_rx_position = true;
if isempty( h_channel(1).rx_position )
    process_rx_position = false;
end

i_bar = 0;                                  % A counter for the progress bar
for i_trk = 1:n_trk                         % Do for each track
    
    % Find the channel-segments in "h_channel" that belong to the current
    % track.
    if only_one_segment
        seg_ind = i_trk;
        n_seg = 1;
    else
        seg_ind = find(strcmp( trk{i_trk} , id_trk ));
        n_seg = numel(seg_ind);
    end
    
    if n_seg > 1
        
        % Check if the track is closed. If it is, the first segment will have an
        % initial position which is > 1.
        if h_channel( seg_ind(1) ).initial_position > 1
            closed = true;
        else
            closed = false;
        end
        
        % Calculate the dimensions of the output channel
        no_tx   = unique( cat( 1 , h_channel( seg_ind ).no_tx) );
        no_rx   = unique( cat( 1 , h_channel( seg_ind ).no_rx) );
        no_path = max( cat( 1 , h_channel( seg_ind ).no_path) ) + 1;
        no_snap = sum( cat( 1 , h_channel( seg_ind ).no_snap) -...
            cat( 1 , h_channel( seg_ind ).initial_position) + 1 );
        
        % Reserve memory for the output data
        coeff = zeros( no_rx , no_tx , no_path , no_snap );
        if h_channel(1).individual_delays
            delay = zeros( no_rx , no_tx , no_path , no_snap );
        else
            delay = zeros( no_path , no_snap  );
        end
        if process_rx_position
            rx_position = zeros( 3,no_snap );
        end
        
        % In order to keep memory requirement of the created data object low,
        % tap positions will be reused. However, we need to store the position
        % of the taps. This is done in "pio".
        pio = zeros(1,no_path);
        pio( 1 : h_channel( seg_ind(1) ).no_path ) = 1 : h_channel( seg_ind(1) ).no_path;
        
        % There might be an phase offset in the LOS component. This will be
        % tracked and compensated.
        phase_offset = ones( no_rx , no_tx );
        
        for i_seg = 1 : n_seg           % Do for each segment
            
            i_bar = i_bar + 1;
            if verbose; m1=ceil(i_bar/n_channel*vb_dots); if m1>m0;
                    for m2=1:m1-m0; fprintf('o'); end; m0=m1; end; end;
            
            % First, we have to determine the indexes of the subsegments.
            % Each segment can be split in three parts: an unused part, an overlapping
            % part and an exclusive part.
            
            % This determines the the start and end-segments for closed tracks
            domerge = true;
            if i_seg < n_seg
                seg_overlap = [ seg_ind(i_seg) , seg_ind(i_seg)+1 ];
            elseif closed
                seg_overlap = [ seg_ind(i_seg) , seg_ind(1) ];
            else
                seg_overlap = seg_ind(i_seg);
                domerge = false;
            end
            
            % === Determine the overlapping area ===
            if domerge
                
                rangem = zeros(3,2);   % The overlapping start / end points
                rangec = zeros(2,2);   % The non-overlapping start / end points
                
                % Determine the end index of the second segment
                rangem(2,2) = h_channel( seg_overlap(2) ).initial_position - 1;
                
                % Determine the start index of the second segment
                rangem(2,1) = floor( rangem(2,2) - (rangem(2,2)-2) * overlap );
                
                % Number of overlapping snapshots
                S = rangem(2,2) - rangem(2,1) + 1;
                
                % Determine the end index of the first segment
                rangem(1,2) = h_channel(seg_overlap(1)).no_snap;
                
                % Determine the start index of the first segment
                rangem(1,1) = rangem(1,2) - S + 1;
                
                % Determine the the non-overlapping area of the first segment
                rangec(1,:) = [ h_channel(seg_overlap(1)).initial_position ,rangem(1,1)-1 ];
                
                % Determine the indexes of the non-overlapping part in the output channel
                if i_seg==1
                    rangec(2,1) = 1;
                else
                    segments_before = seg_ind( seg_ind < seg_ind(i_seg) );
                    rangec(2,1) = sum( cat(1,h_channel( segments_before ).no_snap) -...
                        cat(1,h_channel( segments_before ).initial_position)+1 ) + 1;
                end
                rangec(2,2) = rangec(1,2) - rangec(1,1) + rangec(2,1);
                
                % Determine the indexes of the overlapping part in the output channel
                rangem(3,:) = [rangec(2,2)+1, rangec(2,2)+S];
                
                % Extend the range data structures
                rangem = [rangem(1,1):rangem(1,2) ; rangem(2,1):rangem(2,2) ; rangem(3,1):rangem(3,2) ];
                rangec = [rangec(1,1):rangec(1,2) ; rangec(2,1):rangec(2,2)];
                
            else
                % Here, we only have an exclusive part without overlapping
                rangec = zeros(2,2);   % The non-overlapping start / end points
                
                % Determine the the non-overlapping area of the segment
                rangec(1,:) = [ h_channel(seg_overlap(1)).initial_position ,...
                    h_channel(seg_overlap(1)).no_snap  ];
                segments_before = seg_ind( seg_ind < seg_ind(i_seg));
                
                % Determine the indexes of the non-overlapping part in the output channel
                rangec(2,1) = sum( cat(1,h_channel(segments_before).no_snap) -...
                    cat(1,h_channel(segments_before).initial_position)+1 ) + 1;
                rangec(2,2) = rangec(1,2) - rangec(1,1) + rangec(2,1);
                
                % Extend the range data structure
                rangec = [rangec(1,1):rangec(1,2) ; rangec(2,1):rangec(2,2)];
            end
            % === End ===
            
            
            % Copy the exclusive data from the current segment to the output channel
            coeff( :,:,pio>0,rangec(2,:) ) = h_channel(seg_overlap(1)).coeff( :,:, pio(pio>0) ,...
                rangec(1,:) );
            if h_channel(1).individual_delays
                delay( :,:, pio>0 , rangec(2,:) ) = ...
                    h_channel( seg_overlap(1) ).delay( :,:, pio(pio>0),rangec(1,:) );
            else
                delay( pio>0 , rangec(2,:) ) = ...
                    h_channel( seg_overlap(1) ).delay( pio(pio>0),rangec(1,:) );
            end
            if process_rx_position
                rx_position( :, rangec(2,:) ) = ...
                    h_channel( seg_overlap(1) ).rx_position( :,rangec(1,:) );
            end
            
            % Correct the LOS phase offset of the exclusive part with the value
            % from the previous segment.
            for r = 1:no_rx
                for t = 1:no_tx
                    coeff( r,t,1,rangec(2,:) ) = coeff( r,t,1,rangec(2,:) ) * phase_offset(r,t);
                end
            end
            
            if domerge
                
                % Merge the NLOS components
                L1 = h_channel(seg_overlap(1)).no_path;
                L2 = h_channel(seg_overlap(2)).no_path;
                no_subseg = min([L1,L2]);
                
                % Get the coefficients of the taps
                A = h_channel(seg_overlap(1)).coeff(:,:,1:L1,rangem(1,:)) ;
                B = h_channel(seg_overlap(2)).coeff(:,:,1:L2,rangem(2,:)) ;
                
                
                % Calculate the normalized power of the taps
                pA = reshape( sum(sum(sum(abs(A).^2,1),2),4) ,1,L1);
                pA = pA./sum(pA);
                
                pB = reshape( sum(sum(sum(abs(B).^2,1),2),4) ,1,L2);
                pB = pB./sum(pB);
                
                % Determine the start values for the tap pairing
                [~,oA] = sort(pA(2:end),'descend');
                [~,oB] = sort(pB(2:end),'descend');
                oA = [1 oA+1];
                oB = [1 oB+1];
                
                % "pio" is the position list of the ramped down taps
                % "pin" is the position list of the ramped up taps
                pin = [1 zeros(1,no_path-1)];
                
                % We do not need to merge the subsegments if one of the
                % segments has only a LOS component.
                if no_subseg > 1
                    
                    if optimize
                        % Get the average delays of the taps
                        if h_channel(1).individual_delays
                            dA = h_channel(seg_overlap(1)).delay(:,:,1:L1,rangem(1,:));
                            dA = mean(mean( mean( dA , 1),2),4);
                            dA = reshape(dA,L1,1);
                            
                            dB = h_channel(seg_overlap(2)).delay(:,:,1:L2,rangem(2,:));
                            dB = mean(mean( mean( dB , 1),2),4);
                            dB = reshape(dB,L2,1);
                        else
                            dA = h_channel(seg_overlap(1)).delay(1:L1,rangem(1,:));
                            dA = mean( dA , 2 );
                            
                            dB = h_channel(seg_overlap(2)).delay(1:L2,rangem(2,:));
                            dB = mean( dB , 2 );
                        end
                        
                        % Calculate the initial DS for segment A and segment B
                        dsA = sqrt( pA*dA.^2 - (pA*dA).^2 );
                        dsB = sqrt( pB*dB.^2 - (pB*dB).^2 );
                        
                        pP = [pA,pB];
                        dD = [dA;dB];
                        
                        tmp = no_subseg - 1;
                        ind = (1/(2*tmp) : 1/tmp : 1).';
                        
                        % We want this DS for each subsegment!
                        ds_target = dsA + ( dsB-dsA ) * ind;
                        
                        % the power ramp for LOS and tpas without partner
                        ramp = 0.5*(1+sin((ind-0.5)*pi));
                        
                        % Correct the tap paring iteratively
                        % This is done by first removing a tap from the set and then
                        % changing the ramping order. The position with minimum cost is
                        % saved. This is repeated twice.
                        
                        for cc = 1:2
                            % Find optimal tap ordering for segment A
                            for ca = 2:L1
                                tmp = oA( oA ~= ca );        	% Remove tap from set
                                cost = ones( L1,1 )*Inf;       	% Initialize costs
                                for cb = 2:L1
                                    % Place tap at new position and calculate costs
                                    oAc = [ tmp(1:cb-1) , ca , tmp(cb:end) ];
                                    cost(cb) = merging_cost_fcn( oAc , oB , pP , dD , ds_target , ramp );
                                end
                                [~,cb] = min(cost);             % Look for minimum costs
                                oA = [ tmp(1:cb-1) , ca , tmp(cb:end) ];  % Reorder
                            end
                            
                            % Find optimal tap ordering for segment B
                            for ca = 2:L2
                                tmp = oB( oB ~= ca );        	% Remove tap from set
                                cost = ones( L2,1 )*Inf;       	% Initialize costs
                                for cb = 2:L2
                                    % Place tap at new position and calculate costs
                                    oBc = [ tmp(1:cb-1) , ca , tmp(cb:end) ];
                                    cost(cb) = merging_cost_fcn( oA , oBc , pP , dD , ds_target , ramp );
                                end
                                [~,cb] = min(cost);             % Look for minimum costs
                                oB = [ tmp(1:cb-1) , ca , tmp(cb:end) ];  % Reorder
                            end
                        end
                    end
                    
                    % Calculate the ramp for each sub-interval
                    ramp_length = floor( S/(no_subseg-1) );
                    xl = (1:ramp_length)/(ramp_length+1);
                    xs = 0.5*(1+sin((xl-0.5)*pi));
                    R = [];
                    R(1,1,1,:) = xs;
                    R = repmat(R,no_rx,no_tx);
                    
                    for l = 2 : no_subseg
                        % Calculate the segments
                        ind1 = 1:(l-2)*ramp_length;                       % Before current tap
                        ind2 = (l-2)*ramp_length+1 : (l-1)*ramp_length;   % Current tap
                        ind3 = (l-1)*ramp_length+1 : S;                   % Next taps
                        
                        % Calculate the indices
                        indo = find( pio == oA(l) );               	% Index of the ramp-down tap
                        indn = find( pio == -3 , 1 );               % Index of the ramp-up tap
                        if isempty( indn )
                            indn = find( pio == 0 , 1 );
                        end
                        
                        % Merge the coefficients
                        coeff(:,:,indo,rangem(3,ind1)) = A(:,:,oA(l),ind1);
                        coeff(:,:,indo,rangem(3,ind2)) = A(:,:,oA(l),ind2) .* (1-R);
                        coeff(:,:,indn,rangem(3,ind2)) = B(:,:,oB(l),ind2) .* R;
                        coeff(:,:,indn,rangem(3,ind3)) = B(:,:,oB(l),ind3);
                        
                        % Merge the delays
                        if h_channel(1).individual_delays
                            delay( :,:,indo,rangem(3,ind1) ) = h_channel(seg_overlap(1)).delay( :,:,oA(l),rangem(1,ind1) );
                            delay( :,:,indo,rangem(3,ind2) ) = h_channel(seg_overlap(1)).delay( :,:,oA(l),rangem(1,ind2) );
                            delay( :,:,indn,rangem(3,ind2) ) = h_channel(seg_overlap(2)).delay( :,:,oB(l),rangem(2,ind2) );
                            delay( :,:,indn,rangem(3,ind3) ) = h_channel(seg_overlap(2)).delay( :,:,oB(l),rangem(2,ind3) );
                        else
                            delay( indo,rangem(3,ind1) ) = h_channel(seg_overlap(1)).delay( oA(l),rangem(1,ind1) );
                            delay( indo,rangem(3,ind2) ) = h_channel(seg_overlap(1)).delay( oA(l),rangem(1,ind2) );
                            delay( indn,rangem(3,ind2) ) = h_channel(seg_overlap(2)).delay( oB(l),rangem(2,ind2) );
                            delay( indn,rangem(3,ind3) ) = h_channel(seg_overlap(2)).delay( oB(l),rangem(2,ind3) );
                        end
                        
                        % Update the index list
                        pin( indn ) = oB(l);
                        pio( indo ) = -3;
                        pio( indn ) = -1;
                        
                        % Alternatively, one could use the following code. This
                        % includes a gap between each two taps. However, it
                        % requires to change the number of paths to
                        %    max( cat(1,h_channel(seg_ind).no_path) )+2;
                        % in line 140.
                        
                        % pin( indn ) = oB(l);
                        % pio( pio == -2 ) = -3;
                        % pio( indo ) = -2;
                        % pio( indn ) = -1;
                    end
                end % if no_subseg > 1
                
                % Calculate the ramp for the remaining taps and the LOS
                R = [];                                         % Calculate the Ramp
                x = (1:S)/(S+1);
                R(1,1,1,:) = 0.5*(1+sin((x-0.5)*pi));
                R = repmat(R,no_rx,no_tx);                      % Use same ramp for all links
                
                if L2>L1
                    % Ramp up all remaining taps
                    for l = no_subseg+1 : L2
                        indn = find( pio == 0 , 1 );            % Index of the ramp-up taps
                        
                        coeff(:,:,indn,rangem(3,:)) = B(:,:,oB(l),:) .* R;
                        if h_channel(1).individual_delays
                            delay( :,:,indn,rangem(3,:) ) = h_channel(seg_overlap(2)).delay( :,:,oB(l),rangem(2,:) );
                        else
                            delay( indn,rangem(3,:) ) = h_channel(seg_overlap(2)).delay( oB(l),rangem(2,:) );
                        end
                        
                        pin( indn ) = oB(l);
                        pio( indn ) = -1;
                    end
                else
                    % Ramp down all remaining taps from the old segment
                    for l = no_subseg+1 : L1
                        indo = find( pio == oA(l) );           % Index of the ramp-down tap
                        
                        coeff(:,:,indo,rangem(3,:)) = A(:,:,oA(l),:) .* (1-R);
                        if h_channel(1).individual_delays
                            delay( :,:,indo,rangem(3,:) ) = h_channel(seg_overlap(1)).delay( :,:,oA(l),rangem(1,:) );
                        else
                            delay( indo,rangem(3,:) ) = h_channel(seg_overlap(1)).delay( oA(l),rangem(1,:) );
                        end
                    end
                end
                
                % Merge the LOS components
                A = h_channel(seg_overlap(1)).coeff(:,:,1,rangem(1,:));        % Coefficients of segment 1
                B = h_channel(seg_overlap(2)).coeff(:,:,1,rangem(2,:));        % Coefficients of segment 2
                
                % Correct the phase offset of segment A with the old value
                for r = 1:no_rx
                    for t = 1:no_tx
                        A(r,t,1,:) = A(r,t,1,:) * phase_offset(r,t);
                    end
                end
                
                % Determine new LOS phase offset
                phase_offset = exp( 1j * angle(A(:,:,1,1)) ) ./ exp( 1j * angle(B(:,:,1,1)) );
                
                % Correct the LOS phase offset of the segment B with the new
                % value.
                for r = 1:no_rx
                    for t = 1:no_tx
                        B(r,t,1,:) = B(r,t,1,:) * phase_offset(r,t);
                    end
                end
                
                % Added the coefficients
                coeff( :,:,1,rangem(3,:) ) = A.*(1-R) + B.*(R);
                
                % The LOS delay should be identical for the overlapping area.
                if h_channel(1).individual_delays
                    delay( :,:,1,rangem(3,:) ) = h_channel(seg_overlap(1)).delay( :,:,1,rangem(1,:) );
                else
                    delay( 1,rangem(3,:) ) = h_channel(seg_overlap(1)).delay( 1,rangem(1,:) );
                end
                
                % Process the rx_position for the overlapping part
                if process_rx_position
                    A = h_channel(seg_overlap(1)).rx_position(:,rangem(1,:));        % Positions of segment 1
                    B = h_channel(seg_overlap(2)).rx_position(:,rangem(2,:));        % Positions of segment 2
                    rx_position( :, rangem(3,:) ) = 0.5*( A+B );
                end
                
                % Save the new tap position list for the next segment
                pio = pin;
            end
        end
        
        % If there is only one segment to be merged, an extra tap is generated,
        % but not used. Here, we discard of it.
        if all( reshape( coeff(:,:,end,:) , [] , 1 ) == 0 )
            coeff = coeff(:,:,1:end-1,:);
            if h_channel(1).individual_delays
                delay = delay(:,:,1:end-1,:);
            else
                delay = delay( 1:end-1,: );
            end
        end
        
        % Create an output channel array
        c(i_trk) = channel( coeff,delay );
        c(i_trk).name = trk{i_trk};
        c(i_trk).version = qv;
        if process_rx_position
            c(i_trk).rx_position = rx_position;
        end
        if ~isempty( h_channel(seg_ind(1)).tx_position )
            c(i_trk).tx_position = mean( cat( 2,h_channel(seg_ind).tx_position ),2 );
        end
        
    else % We have only one segment in the current track
        
        % Update progress bar
        i_bar = i_bar + 1;
        if verbose; m1=ceil(i_bar/n_channel*vb_dots); if m1>m0;
                for m2=1:m1-m0; fprintf('o'); end; m0=m1; end; end;
        
        % Copy the original channel
        c(i_trk) = h_channel( seg_ind(1) ).copy;
        c(i_trk).name = trk{i_trk};
    end
end

% Sort channel by name
n_channel = numel(c);
names = {};
for i_channel = 1:n_channel
    names{i_channel} = c(i_channel).name;
end
[~,ind] = sort(names);
c = c(ind);

if verbose
    fprintf('] %5.0f seconds\n',round( etime(clock, tStart) ));
end

end

