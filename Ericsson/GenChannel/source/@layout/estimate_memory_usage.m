function mem = estimate_memory_usage( h_layout, verbose )
%ESTIMATE_MEMORY_USAGE Estimates the required memory usage for the channel generation
%
% QuaDRiGa Copyright (C) 2011-2015 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
%
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

if ~exist( 'verbose' , 'var' )
    verbose = true;
end

% Initialize parameter set structure
% This parses the parameter files, but does not calculate the maps
h_parset = h_layout.create_parameter_sets( -1,1 );

%% Estimate the required Memory for the Maps
mem_maps = 0;
for i_par = 1:numel( h_parset )
    
    % 7 maps plus one in temporary memory for filtering = 8 maps
    % Double precision has 8 bytes per pixel
    tmp = prod( h_parset(i_par).map_size )*8*8;
    
    mem_maps = mem_maps + tmp;
end

if verbose
    disp(' ');
    disp(['Parameter maps:                ',num2str(mem_maps/1024^2 ,'%02.1f' ), ' MB'])
end

%% Memory usage of the channel builder
mem_chan = 0;
mem_cb   = 0;
n_clusters_max = 0;
for i_par = 1:numel( h_parset )
    
    n_tx        = h_parset(i_par).tx_array.no_elements;
    n_clusters  = h_parset(i_par).scenpar.NumClusters;
    n_tx_o      = size( h_parset(i_par).tx_array.coupling,2 );
    
    % Determine size of the tx_pattern
    no_bits = 0;
    if isreal( h_parset(i_par).tx_array.field_pattern_vertical )
        no_bits = no_bits + 8;
    else
        no_bits = no_bits + 16;
    end
    if isreal( h_parset(i_par).tx_array.field_pattern_horizontal )
        no_bits = no_bits + 8;
    else
        no_bits = no_bits + 16;
    end
    tx_pattern_size = h_parset(i_par).tx_array.no_az *...
        h_parset(i_par).tx_array.no_el * n_tx * no_bits;
    
    for i_rx = 1:h_parset(i_par).no_positions
        
        % Determine size of the rx_pattern
        n_rx   = h_parset(i_par).rx_array(i_rx).no_elements;
        n_rx_o = size( h_parset(i_par).rx_array(i_rx).coupling,2 );
        if i_rx == 1 || isequal( h_parset(i_par).rx_array(i_rx) ,...
                h_parset(i_par).rx_array(i_rx-1)  )
            
            no_bits = 0;
            if isreal( h_parset(i_par).rx_array(i_rx).field_pattern_vertical )
                no_bits = no_bits + 8;
            else
                no_bits = no_bits + 16;
            end
            if isreal( h_parset(i_par).rx_array(i_rx).field_pattern_horizontal )
                no_bits = no_bits + 8;
            else
                no_bits = no_bits + 16;
            end
            rx_pattern_size = h_parset(i_par).rx_array(i_rx).no_az *...
                h_parset(i_par).rx_array(i_rx).no_el * n_rx * no_bits;
        end
        
        n_snap = h_parset(i_par).rx_track(i_rx).no_snapshots;
        
        % The size of the output channels for each segment
        tmp_out  = n_rx_o * n_tx_o * n_clusters * n_snap * 8 * 2;
        mem_chan = mem_chan + tmp_out;
        
        % Coefficients for one segment (are stored twice, complex)
        tmp_int1_1  = n_rx * n_tx * n_clusters * n_snap * 16 ;
        tmp_int1_2  = n_rx_o * n_tx_o * n_clusters * n_snap * 16 ;
        
        % Contain subpaths and are stored 3 times
        tmp_int2    = n_rx * n_tx * n_clusters * 20 * 16 * 3;
        
        mem_cb_tmp = tx_pattern_size + rx_pattern_size + tmp_int2 +...
            2*tmp_int1_1 + max([tmp_int1_1,tmp_int1_2]);
        
        % Add space for individual delays
        if h_layout.simpar.drifting_precision > 1
            mem_cb_tmp = mem_cb_tmp + tmp_int1_1/2 + max([tmp_int1_1,tmp_int1_2])/2;
            mem_chan = mem_chan + tmp_out/2;
        end
        
        if mem_cb_tmp > mem_cb
            mem_cb = mem_cb_tmp;
        end
        
        if n_clusters > n_clusters_max
            n_clusters_max = n_clusters + 1;
        end
    end
end

if verbose
    disp(['Internal memory for CB:        ',num2str(mem_cb/1024^2 ,'% 04.1f' ), ' MB']);
    disp(['Output channel file (raw):     ',num2str(mem_chan/1024^2 ,'% 04.1f' ), ' MB']);
end


%% Estimate Memory for the output files
mem_chan_merged = 0;
mem_merger = 0;
for i_tx = 1:h_layout.no_tx
    n_tx        = h_layout.tx_array(i_tx).no_elements;
    n_tx_o      = size( h_layout.tx_array(i_tx).coupling,2 );

    for i_rx = 1:h_layout.no_rx
        n_snap = h_layout.track(i_rx).no_snapshots;
        n_rx_o = size( h_layout.rx_array(i_rx).coupling,2 );
        
        tmp_out =  n_rx_o * n_tx_o * n_clusters_max * n_snap * 8 * 2;
        
        if h_layout.simpar.drifting_precision > 1
            mem_chan_merged = mem_chan_merged + tmp_out*1.5;
        else
            mem_chan_merged = mem_chan_merged + tmp_out;
        end
        
        if tmp_out > mem_merger
            mem_merger = tmp_out;
        end
    end
end


mem1 = mem_maps + mem_cb + mem_chan;
mem2 = mem_maps + mem_chan + mem_merger + mem_chan_merged;

mem = max([mem1 , mem2]);

if verbose
    disp(['Internal memory for  merger:   ',num2str(mem_merger/1024^2 ,'% 04.1f' ), ' MB']);
    disp(['Output channel file (merged):  ',num2str(mem_chan_merged/1024^2 ,'% 04.1f' ), ' MB']);
    disp(' ');
    disp(['Total estimated memory usage:  ',num2str(mem/1024^2 ,'% 04.1f' ), ' MB']);
end

end

