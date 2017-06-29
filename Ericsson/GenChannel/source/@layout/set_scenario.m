function indoor_rx = set_scenario( h_layout, scenario, rx, tx, indoor_frc )
%SET_SCENARIO Assigns scenarios to tracks and segments
%
% This function can be used to assign scenarios to tracks and segments on
% tracks. This takes the distance-dependent LOS probability into account
% for some specific scenarios. Currently, distance-dependent scenario
% selection is available for:
%
%  - 3GPP_3D_UMi
%  - 3GPP_3D_UMa
%  - mmMAGIC_initial_UMi
%  - mmMAGIC_initial_Indoor
%
% Alternatively, you can use all scenarios specified in
% "parameter_set.supported_scenarios".
%
%
% QuaDRiGa Copyright (C) 2011-2016 Fraunhofer Heinrich Hertz Institute
% e-mail: quadriga@hhi.fraunhofer.de
%
% QuaDRiGa is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

no_rx = h_layout.no_rx;
no_tx = h_layout.no_tx;

if ~exist( 'scenario' , 'var' )
    scenario = [];
elseif iscell( scenario ) || ~ischar( scenario )
    error('Scenario must be a string.')
end

if ~exist( 'rx' , 'var' ) || isempty(rx)
    rx = 1 : no_rx;
end

if ~exist( 'tx' , 'var' ) || isempty(tx)
    tx = 1 : no_tx;
end

if ~exist( 'indoor_frc' , 'var' ) || isempty(indoor_frc)
    indoor_frc = 0;
end

% Determine if the user is indoor
tmp = rand( 1, numel(rx) );
indoor_rx = tmp < indoor_frc;

if any( strcmpi( scenario, parameter_set.supported_scenarios ))
    % Set single scenarios
    h_layout.track.set_scenario( scenario );
    
else
    % Determine LOS probabilities and set the scenario accordingly
    for i_rx = 1 : numel(rx)
        rx_ind = rx( i_rx );
        
        % Calculate the Tx-Rx distance for each segment on the track
        segment_index = h_layout.track( rx_ind ).segment_index;
        
        rx_pos_2d = h_layout.track( rx_ind ).positions( 1,segment_index ) +...
            h_layout.track( rx_ind ).initial_position(1) +...
            1j*h_layout.track( rx_ind ).positions( 2,segment_index ) + ...
            1j*h_layout.track( rx_ind ).initial_position(2);
        
        tx_pos_2d = ( h_layout.tx_position( 1,: ) + 1j*h_layout.tx_position( 2,: ) ).';
        
        dist_2d = abs( rx_pos_2d( ones(1,no_tx),: ) -...
            tx_pos_2d( :,ones(1,numel(segment_index)) ) );
        
        % Determine LOS probability
        switch scenario
            case '3GPP_3D_UMi'
                % See: 3GPP TR 36.873 V12.1.0 (2015-03)
                scen = { '3GPP_3D_UMi_LOS', '3GPP_3D_UMi_NLOS',...
                    '3GPP_3D_UMi_LOS_O2I', '3GPP_3D_UMi_NLOS_O2I' };
                
                % Determine the LOS probability for each BS-MT
                p_LOS = min( 18./dist_2d , 1 ) .* (1-exp(-dist_2d/36)) + exp(-dist_2d/36);
                i_LOS = ( rand( size(p_LOS) ) >= p_LOS ) + 1;
                
            case '3GPP_3D_UMa'
                % See: 3GPP TR 36.873 V12.1.0 (2015-03)
                scen = { '3GPP_3D_UMa_LOS', '3GPP_3D_UMa_NLOS',...
                    '3GPP_3D_UMa_LOS_O2I', '3GPP_3D_UMa_NLOS_O2I' };
                
                % Include height-dependency of the user terminals
                h_UT = h_layout.track( rx_ind ).positions( 3,segment_index ) +...
                    h_layout.track( rx_ind ).initial_position(3);
                
                C = zeros( size( dist_2d ));
                
                % Exclude outdoor users from height-dependency
                if indoor_rx(i_rx)
                    g = C;
                    
                    ii = dist_2d > 18;
                    g(ii) = 1.25e-6 .* dist_2d(ii).^2 .* exp( -dist_2d(ii)/150 );
                    
                    ii = h_UT > 13 & h_UT < 23;
                    if any(ii)
                        C( :,ii ) = ones(no_tx,1) * ((h_UT(ii)-13)/10).^1.5;
                    end
                    
                    C( :,h_UT >= 23 ) = 1;
                    C = C .* g;
                end
                
                % Determine the LOS probability for each BS-MT
                p_LOS = ( min( 18./dist_2d , 1 ) .* (1-exp(-dist_2d/63)) + exp(-dist_2d/63) )...
                    .* (1+C);
                i_LOS = ( rand( size(p_LOS) ) >= p_LOS ) + 1;
                
            case 'mmMAGIC_initial_UMi'
                % See: mmMAGIC Deliverable D2.1
                % ICT-671650-mmMAGIC/D2.1
                
                scen = { 'mmMAGIC_initial_UMi_10-80_LOS', 'mmMAGIC_initial_UMi_10-80_NLOS',...
                    'mmMAGIC_initial_UMi_10-80_LOS_O2I', 'mmMAGIC_initial_UMi_10-80_NLOS_O2I' };
                
                % Determine the LOS probability for each BS-MT
                p_LOS = min( 18./dist_2d , 1 ) .* (1-exp(-dist_2d/36)) + exp(-dist_2d/36);
                i_LOS = ( rand( size(p_LOS) ) >= p_LOS ) + 1;
                
            case 'mmMAGIC_initial_Indoor'
                % See: mmMAGIC Deliverable D2.1
                % ICT-671650-mmMAGIC/D2.1
                
                p_LOS = ones(size(dist_2d));
                
                ii = dist_2d > 1.2 & dist_2d < 6.5;
                p_LOS(ii) = exp(-(dist_2d(ii)-1.2)/4.7);
                
                ii = dist_2d >= 6.5;
                p_LOS(ii) = exp(-(dist_2d(ii)-6.5)/32.6)*0.32;
                
                i_LOS = ( rand( size(p_LOS) ) >= p_LOS ) + 1;
                
            otherwise
                error('Scenario is not supported.')
        end
        
        % Set the indoor scenarios
        if indoor_rx(i_rx)
            i_LOS = i_LOS + 2;
        end
        
        tmp = size( h_layout.track( rx_ind ).scenario );
        if tmp(1) == no_tx && tmp(2) == numel(segment_index)
            scen_old = h_layout.track( rx_ind ).scenario;
            
        elseif tmp(1) == 1 && tmp(2) == 1
            scen_old = h_layout.track( rx_ind ).scenario( ones(no_tx,numel(segment_index) ));
            
        elseif tmp(1) == no_tx && tmp(2) == 1
            scen_old = h_layout.track( rx_ind ).scenario( :,ones(1,numel(segment_index) ));
            
        else
            error(['Scenario definition dimension mismatch for Rx ',num2str(rx_ind)]);
        end
        
        scen_current = reshape( scen( i_LOS(:) ) , no_tx,[] );
        
        % Assign scenario to track
        scen_old( tx,: ) = scen_current( tx,: );
        h_layout.track( rx_ind ).scenario = scen_old;
    end
end
end
