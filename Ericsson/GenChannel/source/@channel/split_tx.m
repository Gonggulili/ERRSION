function chan_out = split_tx( h_channel, varargin )
%SPLIT_TX Splits channel arrays based on transmit antenna indices
%   
%   This function can be used to split large transmit antenna arrays into
%   smaller arrays. For example, this can be used to calculate the
%   channels for individual sectors at a BS.
%
%   Example:
%   A channel array has channels from three base stations (BSs). The first
%   and second BS have two sectors, each with two antennas. However, the
%   sector antennas are merged into one array. The third BS has only one
%   sector. To split the channels into five sectors, the following command
%   can be used:    
% 
%       cs = c.split( {1:2,3:4}, {1:2,3:4}, {1:2} );
%
%   Notes: 
%   - The method parses the name-string of the channel objects
%     {channel.name} in order to determine the Tx-Rx relationship. There
%     are two allowed formats: (a) "tx_rx" and (b) "scenario_tx_rx"  
%   - The order of the inputs must match the transmitters in alphabetical
%     order, i.e. the first input corresponds to "Tx01", the second to "Tx02"
%     and so on. This is independent of the order in "layout.tx_name", which
%     might have a different order.   
%   - If only one cell is given as input, but there are several Txs in the
%     channel array, the same sectorization is applied to each one of them.
%   - Outputs are sorted alphabetically according to "tx_rx" (scenario
%     names are ignored). 
%   - If the input array is shaped as [ Rx, Tx ], the output will be shaped
%     as [ Rx, Tx * Sec ]
%
%
% QuaDRiGa Copyright (C) 2011-2015 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
%
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

splt = varargin;

% The number of channel objects
no_c = numel( h_channel );

% Parse names of the channel objects
id_scen = cell(no_c,1);
id_tx = cell(no_c,1);
id_rx = cell(no_c,1);
id_trk = cell(no_c,1);
for n = 1 : no_c
    name = h_channel(n).name;
    p = regexp( name , '_' );
    if numel(p) == 1
        % Format: Tx_Rx
        id_scen{n} = '';
        id_tx{n} = name(1:p(1)-1);
        id_rx{n} = name(p(1)+1:end);
        id_trk{n} = name;
    elseif numel(p) == 2
        % Format: Scenario_Tx_Rx
        id_scen{n} = name(1:p(1)-1);
        id_tx{n} = name(p(1)+1:p(2)-1);
        id_rx{n} = name(p(2)+1:end);
        id_trk{n} = name;
    else
        error('??? Could not parse the channel name string.');
    end
end

% Get the unique tx_names
% The following code implements "unique(id_trk)" but does not sort the
% values in id_trk.

tmp = id_tx;
pt = 1;
while pt < numel(tmp)
    ii = find( strcmp( tmp(pt) , tmp(pt+1:end) ) )+pt;
    tmp = tmp( setdiff( 1:numel(tmp), ii ) );
    pt  = pt + 1;
end
tx_name = tmp;
no_bs = numel( tx_name );
tx_name = sort( tx_name );

% Get the unique rx_names

tmp = id_rx;
pt = 1;
while pt < numel(tmp)
    ii = find( strcmp( tmp(pt) , tmp(pt+1:end) ) )+pt;
    tmp = tmp( setdiff( 1:numel(tmp), ii ) );
    pt  = pt + 1;
end
rx_name = tmp;
no_rx = numel( rx_name );

% Test if input array was sorted

if all( size( h_channel ) == [ no_rx, no_bs ] )
    input_is_sorted = true;
else
    input_is_sorted = false;
end

% Parse input data
if numel( splt ) ~= no_bs
   if  numel( splt ) == 1
       if ~iscell( splt{1} )
           error('??? Inputs must be a cell array.');
       else
           % Copy the data
           for n = 2:no_bs
               splt{n} = splt{1};
           end
       end
   else
        error('??? Number of inputs does not match the number of BSs in the channel array.');
   end
end

for n = 1:numel( splt )
    for m = 1:numel( splt{n} )
        if size( splt{n}{m},1 ) ~= 1
            splt{n}{m} = splt{n}{m}(:)';
        end
        if any( rem( splt{n}{m} , 1 ) ~= 0 )
            error('??? Inputs must be vectors of integer numbers.');
        end
    end
end

chan_out = channel.empty;
sec_cnt  = 1;
for i_bs = 1 : no_bs
    no_sec = numel( splt{i_bs} );
    
    % Do for each input channel
    for i_c = 1 : no_c

        % Check if it matches the input data structure
        if strcmp( id_tx( i_c ), tx_name( i_bs ) )
            
            individual_delays = h_channel(i_c).individual_delays;
            
            for i_sec = 1 : no_sec
                tx_ind = splt{i_bs}{i_sec};
                if ~isempty( id_scen{i_c} )
                    name = [ id_scen{i_c},'_',id_tx{i_c},'s',sprintf('%d', i_sec),'_',id_rx{i_c} ];
                else
                    name = [ id_tx{i_c},'s',sprintf('%d', i_sec),'_',id_rx{i_c} ];
                end
                
                tmp = h_channel(i_c).copy;
                tmp.name  = name;
                tmp.coeff = h_channel(i_c).coeff( :,tx_ind,:,: );
                if individual_delays
                    tmp.delay = h_channel(i_c).delay( :,tx_ind,:,: );
                end
                
                chan_out(sec_cnt) = tmp;
                sec_cnt = sec_cnt + 1;
            end
        end
    end
end

name = cell( numel( chan_out ),1 );
for n = 1 : numel( chan_out )
    p = regexp( chan_out(n).name , '_' );
    if numel(p) == 2
        name{n} = chan_out(n).name(p(1)+1:end);
    else
        name{n} = chan_out(n).name;
    end
end
[~,ind] = sort( name );
chan_out = chan_out(ind);

if input_is_sorted
   chan_out = reshape( chan_out, no_rx , [] );
end

end
