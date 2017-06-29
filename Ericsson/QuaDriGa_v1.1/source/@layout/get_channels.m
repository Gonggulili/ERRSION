function [ h_channel, h_parset, h_cb ] = get_channels( h_layout, sampling_rate, check_parfiles )
%GET_CHANNELS Calculate the channel coefficients
%
%   [ h_channel,h_parset,h_cb ] = GET_CHANNELS generates the channel coefficients.
%       h_channel   is the channel object
%       h_parset    is the parameter_set object
%       h_cb        is the channel_builder object
%
%   Options:
%       sampling_rate (default: 0.01 = 10 ms)
%       The sampling rate in seconds. This parameter is only used if a speed
%       profile is provided by the track object.
%
%       check_parfiles = 0 (default: 1)
%       Disables the parsing of shortnames and the validity-check for the
%       config-files. This is useful, if you know that the parameters in
%       the files are valid. In this case, this saves some execution time.
%
%
% QuaDRiGa Copyright (C) 2011-2015 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
%
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% Parse input variables
if exist( 'sampling_rate','var' ) && ~isempty( sampling_rate )
    if ~( all(size(sampling_rate) == [1 1]) &&...
            isnumeric(sampling_rate) &&...
            isreal(sampling_rate) && sampling_rate > 0 )
        error('??? Invalid sampling interval. The value must be real and > 0.')
    end
else
    sampling_rate = 10e-3; % 10 ms
end

if ~exist( 'check_parfiles' , 'var' )
    check_parfiles = true;
end

verbose = h_layout.simpar.show_progress_bars;

% Check if tracks fulfill the sampling theoreme
samling_limit = h_layout.simpar.wavelength / 2;
sampling_ok = true;
has_speed_profile = false;
for i_rx = 1 : h_layout.no_rx
    if h_layout.track(i_rx).no_snapshots > 1
        [~,dist] = h_layout.track(i_rx).get_length;
        if any( diff(dist) > samling_limit )
            sampling_ok = false;
        end
        if ~isempty( h_layout.track(i_rx).movement_profile )
            has_speed_profile = true;
        end
    end
end

if ~sampling_ok
    warning('QuaDRiGa:layout:get_channels:sampling_ok',...
        'Sample density in tracks does not fulfill the sampling theoreme.');
    
    if has_speed_profile % Change sample density
        h_layout.simpar.sample_density = 2.5;
        for i_rx = 1 : h_layout.no_rx
            if h_layout.track.no_snapshots > 1
                h_layout.track.interpolate_positions( h_layout.simpar.samples_per_meter );
            end
        end
        warning('QuaDRiGa:layout:get_channels:sampling_ok',...
            'Sample density was adjustet to match the sampling theoreme.');
    end
end

% Create LSPs
if isempty( h_layout.track(1).par )
    [ ~, h_parset ] = h_layout.generate_parameters( 0.5, 2, check_parfiles );
else
    h_parset = h_layout.create_parameter_sets( 1,check_parfiles );
end

[h_channel,h_cb] = h_parset.get_channels;
h_channel = h_channel.merge(0.5,1,verbose);

% Get names
n_channel = numel(h_channel);
names = {};
for i_channel = 1:n_channel
    names{i_channel} = h_channel(i_channel).name;
end

% Set the path-gain
for i_rx = 1 : h_layout.no_rx
    if ~isempty( h_layout.track(i_rx).par ) && ~isempty( h_layout.track(i_rx).par.pg )
        
        for i_tx = 1 : h_layout.no_tx
            name = [ h_layout.tx_name{i_tx},'_',h_layout.rx_name{i_rx} ];
            ci   = (i_tx-1)*h_layout.no_rx + i_rx;
            
            % Search all fields
            if ci > numel( h_channel ) || ...
                    isempty( regexp( name ,  h_channel(ci).name, 'once' ) )
                ci = [];
                for i_channel = 1 : n_channel
                    if ~isempty( regexp( name ,  h_channel(i_channel).name, 'once' ) )
                        ci = i_channel;
                    end
                end
            end
            
            if ~isempty(ci)
                if size( h_layout.track(i_rx).par.pg,1 ) == 1
                    h_channel(ci).par.pg = h_layout.track(i_rx).par.pg(1,:);
                else
                    h_channel(ci).par.pg = h_layout.track(i_rx).par.pg(i_tx,:);
                end
            end
        end
    end
end

% Determine, if channel interpolation is needed
values = h_layout.no_rx;
tmp = { h_layout.track.movement_profile };
need_to_interpolate = false;
i_trk = 1;
while ~need_to_interpolate && i_trk <= values
    if ~isempty( tmp{i_trk} )
        need_to_interpolate = true;
    end
    i_trk = i_trk + 1;
end

if need_to_interpolate
    tStart = clock;
    if verbose; fprintf('Interpolate  ['); end; m0=0;
    
    % Apply speed profile, if provided
    channels_done = false(1,n_channel);
    for i_trk = 1 : h_layout.no_rx
        if verbose; m1=ceil(i_trk/values*50); if m1>m0; for m2=1:m1-m0; fprintf('o'); end; m0=m1; end; end;
        
        trk = h_layout.track(i_trk);
        
        if ~isempty( trk.movement_profile )
            pos_snap = [0,cumsum(abs(diff( trk.positions(1,:) + 1j*trk.positions(2,:) )))];
            dist = trk.interpolate_movement( sampling_rate );
            length = trk.get_length;
            
            for i_channel = 1 : n_channel
                if ~channels_done( i_channel )
                    if ~isempty( regexp( names{i_channel} ,  trk.name, 'once' ) )
                        par = h_channel(i_channel).par;
                        
                        % Interpolate path gain
                        if isfield( par,'pg' ) && size(par.pg,2) == h_channel(i_channel).no_snap
                            par.pg = spline( pos_snap , par.pg , dist );
                        end
                        
                        h_channel(i_channel) = h_channel(i_channel).interpolate( dist, 'spline' );
                        h_channel(i_channel).par = par;
                        channels_done( i_channel ) = true;
                    end
                end
            end
        end
    end
    
    if verbose
        fprintf('] %5.0f seconds\n',round( etime(clock, tStart) ));
    end
end

% Reshape the channel object
if numel(h_channel) == h_layout.no_rx * h_layout.no_tx
    h_channel = reshape( h_channel, h_layout.no_rx , h_layout.no_tx );
end

end
